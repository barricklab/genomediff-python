# -*- coding: utf-8 -*-
"""
 * @Date: 2024-12-27 17:48:21
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-10 21:40:56
 * @FilePath: /pymummer/genomediff/__init__.py
 * @Description:
 modified from https://github.com/biosustain/genomediff-python
"""
# """
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
        for i, record in parse(fsock):
            if isinstance(record, Metadata):
                metadata.set(record.name, record.value)
            elif isinstance(record, str):
                comments[i] = record
            else:
                records.set(record)
        if not hasattr(gdfile, "read"):
            fsock.close()
        return cls(metadata, records, comments)

    def write(self, fsock):
        raise NotImplementedError()

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

    def __getitem__(self, key):
        return self.records[key]
