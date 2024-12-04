import random
import pulp

from typing import Sequence

from rna3db.utils import PathLike, read_json, write_json


def find_optimal_components(
    components: Sequence[int], bins: Sequence[int], verbose: bool = False
) -> Sequence[set[int]]:
    """Function used to find optimal placement of components into
    training/testing sets.

    We use an ILP formulation that is very similar to the classic ILP
    formulation of the bin packing problem.

    Args:
        components (Sequence[int]): list of component sizes
        bins (Sequence[int]): list of bin sizes
        verbose (bool): whether to print verbose output
    Returns:
        Sequence[set[int]]: list of sets, where each set contains the indices
            of the components that go into that bin
    """

    n, k = len(components), len(bins)

    # set up problem
    p = pulp.LpProblem("OptimalComponentSolver", pulp.LpMinimize)
    x = pulp.LpVariable.dicts(
        "x", ((i, j) for i in range(n) for j in range(k)), cat="Binary"
    )
    deviation = pulp.LpVariable.dicts(
        "d", (j for j in range(k)), lowBound=0, cat="Continuous"
    )

    # we want to minimise total "deviation"
    # (deviation is the total sum of the difference between target bins and found bins)
    p += pulp.lpSum(deviation[j] for j in range(k))

    # components can go into exactly one bin
    for i in range(n):
        p += pulp.lpSum(x[(i, j)] for j in range(k)) == 1, f"AssignComponent_{i}"

    # deviation constraints (to handle abs)
    for j in range(k):
        total_weight_in_bin = pulp.lpSum(components[i] * x[(i, j)] for i in range(n))
        p += total_weight_in_bin - bins[j] <= deviation[j], f"DeviationPos_{j}"
        p += bins[j] - total_weight_in_bin <= deviation[j], f"DeviationNeg_{j}"

    # solve ILP problem with PuLP
    p.solve(pulp.PULP_CBC_CMD(msg=int(verbose)))

    # extract solution in sensible format
    sol = [set() for i in range(k)]
    for i in range(k):
        for j in range(n):
            if pulp.value(x[(j, i)]) == 1:
                sol[i].add(j)

    return sol


def split(
    input_path: PathLike,
    output_path: PathLike = None,
    splits: Sequence[float] = [0.7, 0.0, 0.3],
    split_names: Sequence[str] = ["train_set", "valid_set", "test_set"],
    shuffle: bool = False,
    force_zero_last: bool = False,
):
    """A function that splits a JSON of components into a train/test set.

    The split is done by assigning components into the training set until a
    specified training set split percentage (train_size) is met. This is done
    starting with the largest component.

    Args:
        input_path (PathLike): path to JSON containing components
        output_path (PathLike): path to output JSON
    """
    if sum(splits) != 1.0:
        raise ValueError("Sum of splits must equal 1.0.")

    if len(splits) != len(split_names):
        raise ValueError("Number of splits must match number of split names.")

    cluster_json = read_json(input_path)

    # get lengths of the components, and mapping from idx to keys
    keys, lengths = [], []
    for k, v in cluster_json.items():
        if force_zero_last and k == "component_0":
            continue
        keys.append(k)
        lengths.append(len(v))

    # calculate actual bin capacities
    # rounding is probably close enough
    bins = [round(sum(lengths) * ratio) for ratio in splits]

    # create output dict
    output = {k: {} for k in split_names}

    # force `component_0` into the last bin
    if force_zero_last:
        if bins[-1] < len(cluster_json["component_0"]):
            print(
                "ERROR: cannot force `component_0` into the last bin. Increase the last bin size."
            )
            raise ValueError
        bins[-1] -= len(cluster_json["component_0"])
        output[split_names[-1]]["component_0"] = cluster_json["component_0"]
        del cluster_json["component_0"]

    if shuffle:
        L = list(zip(keys, lengths))
        random.shuffle(L)
        keys, lengths = zip(*L)
        keys, lengths = list(keys), list(lengths)

    # find optimal split with ILP
    sol = find_optimal_components(lengths, bins)

    # write output to dict
    for idx, name in enumerate(split_names):
        for k in sorted(sol[idx]):
            k = keys[k]
            output[name][k] = cluster_json[k]

    if output_path:
        write_json(output, output_path)

    return output
