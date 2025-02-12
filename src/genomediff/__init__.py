# -*- coding: utf-8 -*-
from pathlib import Path
from typing import TextIO

from .parser import Metadata, MetadataContainer, load_json_ref_cov, parse
from .records import RecordCollection


class GenomeDiff:
    def __init__(
        self,
        metadata: MetadataContainer,
        records: RecordCollection,
        comments: dict[int, str] | None = None,
    ):
        self._metadata = metadata
        self._records = records
        self._comments = comments if comments else {}

    @property
    def metadata(self):
        return self._metadata

    @property
    def records(self):
        return self._records

    @classmethod
    def read(cls, gdfile: str | Path | TextIO):
        records = RecordCollection.new()
        comments: dict[int, str] = {}
        if hasattr(gdfile, "read"):
            metadata = MetadataContainer(getattr(gdfile, "name", ""))
            fsock: TextIO = gdfile  # type: ignore[assignment]
        else:
            metadata = MetadataContainer(gdfile)
            fsock = open(gdfile, "r")
        self = cls(metadata, records, comments)
        for i, record in parse(fsock):
            if isinstance(record, Metadata):
                metadata.set(record.name, record.value)
            elif isinstance(record, str):
                comments[i] = record
            else:
                record.document = self
                records.set(record)
        if not hasattr(gdfile, "read"):
            fsock.close()
        return self

    def write(self, fsock):
        for l in self.metadata.lines:
            print(l, file=fsock)
        for record in self.records:
            print(str(record), file=fsock)

    @property
    def cov_summary(self):
        if outdir := self.metadata.output:
            return load_json_ref_cov(outdir)
        raise AttributeError("No output directory found in metadata")

    @property
    def mutations(self):
        return list(self.records.mutation)

    @property
    def evidence(self):
        return list(self.records.evidence)

    @property
    def validation(self):
        return list(self.records.validation)

    def __getitem__(self, key):
        return self.records[key]
