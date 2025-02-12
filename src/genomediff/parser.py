import json
import re
from collections import OrderedDict
from pathlib import Path
from time import strptime
from typing import Final, Iterable, NamedTuple
from warnings import warn

import pandas as pd

from .records import RecordEnum

METADATA_PATTERN: Final = re.compile(r"^#=([^\s]+)\s+(.*)$")
MUTATION_PATTERN: Final = re.compile(
    r"^(?P<type>[A-Z]{2,4})"
    r"\t(?P<id>\d+|\.)"
    r"\t((?P<parent_ids>\d+(,\s*\d+)*)|\.?)"
    r"\t(?P<extra>.+)?$"
)


class Metadata(NamedTuple):
    name: str
    value: str

    def __repr__(self):
        return f"#={repr(self.name)} {repr(self.value)}"


class MetadataContainer:
    """https://github.com/barricklab/breseq/wiki/GenomeDiff-File-Format#metadata-lines

    Lines beginning with **#=<name> <value>** are interpreted as metadata. (Thus, the first line is assigning a metadata item named GENOME_DIFF a value of 1.0.) Names cannot include whitespace characters. Values may include whitespace characters. In most cases, the values for lines with the same name are concatenated with single spaces added between them or interpreted as a list.

    Some of these metadata fields are used to name and sort samples by |gdtools| COMPARE, to find relevant files by |gdtools| RUNFILE, and by other utilities.
    """

    def __init__(self, file_path: str | Path):
        self._dict: OrderedDict[str, list[str]] = OrderedDict()
        self.file_path = Path(file_path)

    @property
    def lines(self):
        for k, vs in self._dict.items():
            for v in vs:
                yield Metadata(k, v)

    @property
    def GENOME_DIFF(self):
        """The first line of the file must define that this is a |GD| file and the version of the file specification used::

        #=GENOME_DIFF 1.0
        """
        (version,) = self._dict["GENOME_DIFF"]
        return version

    @property
    def gd_version(self):
        return self.GENOME_DIFF

    def check_multi_values(self, key: str):
        if (values := self._dict.get(key)) is None:
            return None
        if len(values) > 1:
            warn(f"Unexpected multiple values for {key}: {values}")
        return values[-1]

    # Common but optional metadata fields include:
    @property
    def TITLE(self):
        """The name of the sample.
        If this field is not provided, the name of the |GD| file (removing the .gd suffix) is used for this field.
        """
        if "TITLE" not in self._dict:
            return self.file_path.with_suffix("").name
        return self.check_multi_values("TITLE")

    @property
    def AUTHOR(self):
        """Name of person who curated the |GD| file."""
        return self.check_multi_values("AUTHOR")

    @property
    def PROGRAM(self):
        """Name and version of software program that generated the |GD| file."""
        return self.check_multi_values("PROGRAM")

    @property
    def CREATED(self):
        """Date on which the |GD| file was created."""
        t_ = self.check_multi_values("CREATED")
        if t_ is None:
            return None
        return strptime(t_, "%H:%M:%S %d %b %Y")

    @property
    def TIME(self):
        """Time point the sample is from, in days, generations, or any other unit of measurement. Ex: 1, 2, 15000"""
        return self.check_multi_values("TIME")

    @property
    def POPULATION(self):
        """Name/designation for the population the sample is from. Ex: Ara–3 / MA-1"""
        return self.check_multi_values("POPULATION")

    @property
    def TREATMENT(self):
        """Experimental treatment group for this population. Ex: LB medium / LTEE"""
        return self.check_multi_values("TREATMENT")

    @property
    def CLONE(self):
        """Name/designation for a clonal isolate, Ex: A, B, REL10863"""
        return self.check_multi_values("CLONE")

    @property
    def REFSEQ(self):
        """Location of the reference sequence file. Ex: /here/is/an/absolute/path/to/the/file.gb"""
        return [Path(i) for i in self._dict.get("REFSEQ", [])]

    @property
    def ADAPTSEQ(self):
        """Location of the adaptor sequence file. Ex: relative/path/to/the/adaptors.fa"""
        return [Path(i) for i in self._dict.get("ADAPTSEQ", [])]

    @property
    def READSEQ(self):
        """Location of the read sequence file. Ex: https://place.org/url/for/file/download.fastq"""
        return [Path(i) for i in self._dict.get("READSEQ", [])]

    @property
    def COMMAND(self):
        """Command line used to generate the |GD| file."""
        return self.check_multi_values("COMMAND")

    @property
    def output(self):
        if cmd := self.COMMAND:
            cmds = cmd.split()
            index = -1
            if (index := cmds.index("-o")) != -1:
                return Path(cmds[index + 1])
            if (index := cmds.index("--output")) != -1:
                return Path(cmds[index + 1])

    def set(self, name: str, value: str):
        self._dict.setdefault(name, []).append(value)

    @classmethod
    def from_dict(cls, d: dict[str, str | list]):
        self = cls("")
        for k, v in d.items():
            if isinstance(v, list):
                for i in v:
                    self.set(k, str(i))
            else:
                self.set(k, str(v))
        return self

    def __eq__(self, value):
        if isinstance(value, MetadataContainer):
            return dict(self._dict) == dict(value._dict)
        if isinstance(value, dict):
            return self == self.from_dict(value)


def read_line(line: str):
    if line.startswith("#"):
        match = METADATA_PATTERN.match(line)
        if match:
            return Metadata(*match.group(1, 2))
        """
        Lines beginning with whitespace and # are comments. Comments may not occur at the end of a data line.
        """
        return line
    else:
        match = MUTATION_PATTERN.match(line)
        if match:
            record_type = match.group("type")
            assert record_type in RecordEnum._member_names_
            record = RecordEnum[record_type].value(
                match.group("id"),
                match.group("parent_ids"),
                match.group("extra"),
            )
            return record
    return


def parse(fsock: Iterable[str]):
    for line_ptr, line in enumerate(fsock):
        if not line:
            continue
        record = read_line(line)
        if record:
            yield line_ptr, record
        else:
            raise Exception(f"Could not parse line #{line_ptr}: {line}")


def GenomeDiffParser(fsock: Iterable[str]):
    return (j for i, j in parse(fsock=fsock))


def load_json_ref_cov(
    breseq_output_dir: Path | str, relative_file: str = "output/summary.json"
):
    breseq_json = (
        (Path(breseq_output_dir) / relative_file)
        if relative_file
        else breseq_output_dir
    )
    with open(breseq_json) as ji:
        j = json.load(ji)
    df = (
        pd.DataFrame(j["references"]["reference"])
        .T.rename_axis(index="reference")
        .assign(breseq_output_dir=breseq_json)
    )
    return df
