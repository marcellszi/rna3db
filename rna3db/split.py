from rna3db.utils import PathLike, read_json, write_json


def split(
    input_path: PathLike,
    output_path: PathLike,
    train_size: float = 0.7,
    force_zero_test: bool = True,
):
    """A function that splits a JSON of components into a train/test set.

    The split is done by assigning components into the training set until a
    specified training set split percentage (train_size) is met. This is done
    starting with the largest component.

    Args:
        input_path (PathLike): path to JSON containing components
        output_path (PathLike): path to output JSON
        train_size (float): percentage of data to use as training set
        force_zero_test (bool): whether to force component_0 into the test set
    """
    # read json
    cluster_json = read_json(input_path)
    # count number of repr sequences
    total_repr_clusters = sum(len(v) for v in cluster_json.values())

    # figure out which components need to go into training set
    train_components = set()
    train_set_length = 0
    i = 1 if force_zero_test else 0
    while train_set_length / total_repr_clusters < train_size:
        # skip if it's not a real component (should only happen with 0)
        if f"component_{i}" not in cluster_json:
            i += 1
            continue
        train_components.add(f"component_{i}")
        train_set_length += len(cluster_json[f"component_{i}"].keys())
        i += 1

    # test_components are just total-train_components
    test_components = set(cluster_json.keys()) - train_components

    # actually build JSON
    output = {"train_set": {}, "test_set": {}}
    for k in sorted(train_components):
        output["train_set"][k] = cluster_json[k]
    for k in sorted(test_components):
        output["test_set"][k] = cluster_json[k]

    write_json(output, output_path)
