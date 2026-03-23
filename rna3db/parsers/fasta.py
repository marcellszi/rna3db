from pathlib import Path


class FASTA:
    """Parsed FASTA file storing parallel lists of headers and sequences."""

    def __init__(self, headers: list[str], sequences: list[str]):
        """
        Args:
            headers (list[str]): Sequence identifiers (without the ``>`` prefix).
            sequences (list[str]): Corresponding sequences.
        """
        self.headers = list(headers)
        self.sequences = list(sequences)

    def __len__(self) -> int:
        return len(self.headers)

    def __iter__(self):
        return zip(self.headers, self.sequences)

    def __getitem__(self, idx) -> tuple[str, str]:
        return (self.headers[idx], self.sequences[idx])

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FASTA):
            return NotImplemented
        return self.headers == other.headers and self.sequences == other.sequences

    def unpack(self) -> tuple[list[str], list[str]]:
        """Return headers and sequences as a tuple of two lists.

        Example:
            >>> headers, sequences = FASTA.read("seqs.fa").unpack()
        """
        return self.headers, self.sequences

    def __repr__(self) -> str:
        return f"FASTA(n={len(self)})"

    @classmethod
    def read(cls, path: Path, force_gzip: bool = False) -> "FASTA":
        """Parse a FASTA file.

        Supports multi-line sequences.

        Args:
            path (Path): Path to input FASTA file.
            force_gzip (bool, optional): If True, will attempt to read the file as a
                gzip file.

        Returns:
            FASTA: Parsed FASTA file.
        """
        if Path(path).suffix == ".gz" or force_gzip:
            import gzip

            reader = gzip.open(path, "rt")
        else:
            reader = open(path, "r")

        headers = []
        sequences = []
        current_header = None
        current_sequence = ""

        with reader as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        headers.append(current_header)
                        sequences.append(current_sequence)
                    current_header = line[1:]
                    current_sequence = ""
                elif line.startswith("#") or not line:
                    continue
                else:
                    current_sequence += line
            if current_header is not None:
                headers.append(current_header)
                sequences.append(current_sequence)

        return cls(headers, sequences)

    def write(self, path: Path):
        """Write FASTA records to a file.

        Args:
            path (Path): Path to write FASTA file to.
        """
        with open(path, "w") as f:
            for header, sequence in zip(self.headers, self.sequences):
                f.write(f">{header}\n{sequence}\n")
