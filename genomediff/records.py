# -*- coding: utf-8 -*-
"""
 * @Date: 2024-12-27 17:50:21
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-10 21:41:15
 * @FilePath: /pymummer/genomediff/records.py
 * @Description:
"""
# """

import re
from collections import OrderedDict
from enum import Enum
from typing import Final, Iterable, Literal, NamedTuple

from git import TYPE_CHECKING

STRAND_MAPPING_PATTERN: Final = re.compile(r"^(\d+)/(\d+)$")


class ReadEvidenceStrand(NamedTuple):
    pos: int = 0
    neg: int = 0

    def __str__(self):
        return f"{self.pos}/{self.neg}"


def _convert_value(value: str):
    if value == "." or value == "" or value is None:
        return None
    for type_ in (int, float):
        try:
            return type_(value)
        except ValueError:
            pass
    if (match := STRAND_MAPPING_PATTERN.match(value)) is not None:
        return ReadEvidenceStrand(int(match.group(1)), int(match.group(2)))
    return value


class DataMeta(type):
    types: dict[str, str] = {}

    def __init__(cls, name, bases, attrs):
        super().__init__(name, bases, attrs)
        if attrs.get("DataItem", NotImplemented) != NotImplemented:
            cls.types[name] = bases[0].type_class


class DataAbstract(metaclass=DataMeta):
    type_class: str = NotImplemented

    def __init__(self, evidence_id: str, parent_ids: str, extra: str):
        """Data lines describe either a mutation or evidence from an analysis that can potentially support a mutational event. Data fields are tab-delimited. Each line begins with several fields containing information common to all types, continues with a fixed number of type-specific fields, and ends with an arbitrary number of name=value pairs that store optional information.

        1. **type** *<string>*

           type of the entry on this line.

        2. **id or evidence-id** *<uint32>*

           For evidence and validation lines, the id of this item. For mutation lines, the ids of all evidence or validation items that support this mutation. May be set to '.' if a line was manually edited.

        3. **parent-ids** *<uint32>*

           ids of evidence that support this mutation. May be set to '.' or left blank.
        """
        self.evidence_id = _convert_value(evidence_id)
        self.parent_ids = (
            [
                parent_id
                for parent_ids in parent_ids.split(",")
                if (parent_id := _convert_value(parent_ids)) is not None
            ]
            if _convert_value(parent_ids) is not None
            else []
        )
        self.parse_extra(extra)

    DataItem = NotImplemented

    def load_data_item(self, extra: str):
        _fields = self.DataItem._fields
        extra_split = extra.split("\t")
        self.dataitem = self.DataItem(  # type: ignore[misc]
            *(self.DataItem.__annotations__[k](v) for k, v in zip(_fields, extra_split))
        )
        return extra_split[len(_fields) :]

    def parse_extra(self, extra: str):
        extra2dict = self.load_data_item(extra)
        self.extra = OrderedDict(
            (k, _convert_value(v))
            for k, v in (e.split("=", 1) for e in extra2dict if e)
        )

    @property
    def type(self):
        return self.__class__.__name__

    @property
    def id(self):
        """alias"""
        return self.evidence_id

    def __getattr__(self, item):
        if item in self.dataitem._fields:
            return self.dataitem.__getattribute__(item)
        return self.extra[item]

    def __repr__(self):
        return "\t".join(
            [
                self.type,
                str(self.id),
                " ".join([str(i) for i in self.parent_ids]),
            ]
            + [str(i) for i in self.dataitem]
            + [f"{k}={v}" for k, v in self.extra.items()]
        )

    def to_dict(self) -> dict:
        return {"type": self.type, "id": self.id} | self.dataitem._asdict() | self.extra

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.type == other.type
                and self.id == other.id
                and self.dataitem == other.dataitem
                and self.extra == other.extra
            )
        return False

    def __le__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.type == other.type
                and self.id == other.id
                and self.dataitem == other.dataitem
                and self.extra
                == {k: v for k, v in other.extra.items() if k in self.to_dict()}
            )
        return False

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self <= other and self != other


# fmt: off
class DataMutationAbs  (DataAbstract):
    type_class = "mutation"
    def range(self):
        return range(self.dataitem.position, self.dataitem.position + self.dataitem.size)

class DataEvidenceAbs  (DataAbstract): type_class = "evidence"
class DataValidationAbs(DataAbstract): type_class = "validation"
# fmt: on

if TYPE_CHECKING:
    SeqStrand = Literal["+", "-"]
else:

    class SeqStrand(str):
        ARGS = {"1": "+", "-1": "-"}

        def __new__(cls, value):
            if value in cls.ARGS.values():
                return super().__new__(cls, value)
            value = str(value)
            if value in cls.ARGS:
                return super().__new__(cls, cls.ARGS[value])
            raise ValueError(f"Invalid strand value: {value}")


# fmt: off
class SNP(DataMutationAbs):
    """Base substitution mutation
    - seq_id   id of reference sequence fragment containing mutation, evidence, or validation.
    - position position in reference sequence fragment of base to replace.
    - new_seq  new base at position.
    """
    class DataItem(NamedTuple): seq_id: str; position: int; new_seq: str
    dataitem: DataItem
    def range(self): return range(self.dataitem.position, self.dataitem.position + 1)
    # SNP("1", "1", "1\t1\t1").dataitem.new_seq
class SUB(DataMutationAbs):
    """Multiple base substitution mutation
    - seq_id
    - position position in reference sequence of the first base that will be replaced.
    - size     number of bases *after* the specified reference position to replace with **new_seq**.
    - new_seq  new bases to substitute.
    """
    class DataItem(NamedTuple): seq_id: str; position: int; size: int; new_seq: str
    dataitem: DataItem
class DEL(DataMutationAbs):
    """Deletion mutation
    - seq_id
    - position position in reference sequence fragment of first deleted base.
    - size     number of bases deleted in reference.

    * **mediated=** *<mobile_element_family>*
        This deletion appears to be mediated by a molecular event involving a mobile element such as a transposon. A copy of the mobile element is found on the boundary of the deleted region and a new junction at the opposite end of the deletion matches the end of the mobile element.
    * **between=** *<repeat_family>*
        This deletion appears to result from homologous recombination or polymerase slipping between two existing copies of the same genomic repeat (e.g. tRNA, IS element) in the genome. One copy of the repeat is deleted by this event.
    * **repeat_seq=** *<string>*, **repeat_length=** *<uint32>*, **repeat_ref_num=** *<uint32>*, **repeat_new_copies=** *<uint32>*
        This deletion is in a short sequence repeat consisting of tandem copies of **repeat_seq** repeated **repeat_ref_num** times in the ancestor and **repeat_new_copies** after a mutation.  To be annotated in this way the copy of the repeat in the reference genome must consist of at least two repeat copies and have a length of five of more total bases (**repeat_length** × **repeat_ref_num** ≥ 5).
    """
    class DataItem(NamedTuple): seq_id: str; position: int; size: int
    dataitem: DataItem
class INS(DataMutationAbs):
    """Insertion mutation
    - seq_id
    - position position in reference sequence fragment. New bases are inserted *after* this position.
    - new_seq  new bases to be inserted in the reference.

    * **repeat_seq=** *<string>*, **repeat_length=** *<uint32>*, **repeat_ref_num=** *<uint32>*, **repeat_new_copies=** *<uint32>*
        This insertion is in a short sequence repeat consisting of tandem copies of **repeat_seq** repeated **repeat_ref_num** times in the ancestor and **repeat_new_copies** after a mutation.  To be annotated in this way the copy of the repeat in the reference genome must consist of at least two repeat copies and have a length of five of more total bases (**repeat_length** × **repeat_ref_num** ≥ 5).
    * **insert_position=** *<uint32>*
        Used when there are multiple insertion events after the same reference base to order the insertions. This typically happens in polymorphism mode and when manually breaking up an insertion of bases into distinct mutational events when this is supported by phylogenetic information. Numbering of insert positions begins with 1.
    """
    class DataItem(NamedTuple): seq_id: str; position: int;            new_seq: str
    dataitem: DataItem
    def range(self): return range(self.dataitem.position, self.dataitem.position + 1)
class MOB(DataMutationAbs):
    """Mobile element insertion mutation
    - seq_id
    - position         position in reference sequence fragment of the first duplicated base at the target site.
    - repeat_name      name of the mobile element. Should correspond to an annotated **repeat_region** or **mobile_element** feature in the reference sequence.
    - strand           strand of mobile element insertion.
    - duplication_size number of target site bases duplicated during insertion of the mobile element, beginning with the specified reference position. If the value of this field is negative, then it indicates that the absolute value of this number of bases were deleted at the target site beginning with the specified position. If the value of this field is zero, then the there were no duplicated bases, and the mobile element was inserted after the specified base position.

    * **del_start=** *<uint32>*, **del_end=** *<uint32>*
        Delete this many bases from the start or end of the inserted mobile element. This deletion occurs with respect to the top strand of the genome after the element is flipped to the orientation with which it will be inserted.
    * **ins_start=** *<string>*, **ins_end=** *<string>*
        Append the specified bases to the start or end of the inserted mobile element. These insertions occur after any deletions and will be inside of any duplicated target site bases.
    * **mob_region** = *<seq_id:start-end >*
        Use the existing copy of the mobile element specified as a seq_id:start-end region to apply this mutation. Useful when different annotated members of a mobile element family have slightly different sequences.
    """
    class DataItem(NamedTuple): seq_id: str; position: int;                         repeat_name: str; strand: SeqStrand; duplication_size: int
    dataitem: DataItem
    def range(self): raise NotImplementedError
class AMP(DataMutationAbs):
    """Amplification mutation
    - seq_id
    - position         position in reference sequence fragment.
    - size             number of bases duplicated starting with the specified reference position.
    - new_copy_number  new number of copies of specified bases.

    * **between=** *<repeat_family>*
        This amplification appears to result from homologous recombination or polymerase slipping between two existing copies of the same genomic repeat (e.g. tRNA, IS element) in the genome. This repeat appears on the boundary of each copy of the specified region.
    * **mediated=** *<repeat_family>*, *mediated_strand=** *<1/-1>*
        This amplification is mediated by a simultaneous new insertion of a mobile element (or other  repeat element). New copies of the inserted element are added in the specified strand orientation between each new copy of the amplified region. Both of these attributes must be specified for the mutation.
    * **mob_region** = *<seq_id:start-end >*
        Only valid for 'mediated' amplifications. Use the existing copy of the mobile element specified as a seq_id:start-end region to apply this mutation. Useful when different annotated members of a mobile element family have slightly different sequences.
    """
    class DataItem(NamedTuple): seq_id: str; position: int; size: int;              new_copy_number: int
    dataitem: DataItem
class CON(DataMutationAbs):
    """Gene conversion mutation
    - seq_id
    - position position in reference sequence fragment that was the target of gene conversion from another genomic location.
    - size     number of bases to replace in the reference genome beginning at the specified position.
    - region   Region in the reference genome to use as a replacement.
    """
    class DataItem(NamedTuple): seq_id: str; position: int; size: int;              region: str
    dataitem: DataItem
class INV(DataMutationAbs):
    """Inversion mutation
    - seq_id
    - position position in reference sequence fragment.
    - size     number of bases in inverted region beginning at the specified reference position.
    """
    class DataItem(NamedTuple): seq_id: str; position: int; size: int
    dataitem: DataItem


class RA(DataEvidenceAbs):
    """Read alignment evidence
    - seq_id
    - position        position in reference sequence fragment.
    - insert_position number of bases inserted after the reference position to get to this base. An value of zero refers to the base. A value of 5 means that this evidence if for the fifth newly inserted column after the reference position.
    - ref_base        base in the reference genome.
    - new_base        new base supported by read alignment evidence.
    """
    class DataItem(NamedTuple): seq_id: str; position: int;                         insert_position: int; ref_base: str; new_base: str
    dataitem: DataItem
class MC(DataEvidenceAbs):
    """Missing coverage evidence
    - seq_id
    - start       start position in reference sequence fragment.
    - end         end position in reference sequence of region.
    - start_range number of bases to offset *after* the **start position** to define the upper limit of the range where the start of a deletion could be.
    - end_range   number of bases to offset *before* the **end position** to define the lower limit of the range where the start of a deletion could be.
    """
    class DataItem(NamedTuple): seq_id: str; start: int; end: int;                  start_range: int; end_range: int
    dataitem: DataItem
class JC(DataEvidenceAbs):
    """New junction evidence
    - side_1_seq_id   id of reference sequence fragment containing side 1 of the junction.
    - side_1_position position of side 1 at the junction boundary.
    - side_1_strand   direction that side 1 continues matching the reference sequence
    - side_2_seq_id   id of reference sequence fragment containing side 2 of the junction.
    - side_2_position position of side 2 at the junction boundary.
    - side_2_strand   direction that side 2 continues matching the reference sequence.
    - overlap         Number of bases that the two sides of the new junction have in common.
    """
    class DataItem(NamedTuple): side_1_seq_id: str; side_1_position: int; side_1_strand: SeqStrand; side_2_seq_id: str; side_2_position: int; side_2_strand: SeqStrand; overlap: int
    dataitem: DataItem
class CN(DataEvidenceAbs):
    """Copy number variation evidence
    - seq_id
    - start       start position in reference sequence fragment.
    - end         end position in reference sequence of region.
    - copy_number number of copies of the region in the reference genome.
    """
    class DataItem(NamedTuple): seq_id: str; start: int; end: int;                  copy_number: int
    dataitem: DataItem
class UN(DataEvidenceAbs):
    """Unknown base evidence
    - seq_id
    - start  start position in reference sequence of region.
    - end    end position in reference sequence of region.
    """
    class DataItem(NamedTuple): seq_id: str; start: int; end: int
    dataitem: DataItem


class CURA(DataValidationAbs):
    class DataItem(NamedTuple): expert: str
    dataitem: DataItem
class FPOS(DataValidationAbs):
    class DataItem(NamedTuple): expert: str
    dataitem: DataItem
class PHYL(DataValidationAbs):
    class DataItem(NamedTuple): gd: str
    dataitem: DataItem
class TSEQ(DataValidationAbs):
    class DataItem(NamedTuple): seq_id: str; primer1_start: int; primer1_end: int; primer2_start: int; primer2_end: int
    dataitem: DataItem
class PFLP(DataValidationAbs):
    class DataItem(NamedTuple): seq_id: str; primer1_start: int; primer1_end: int; primer2_start: int; primer2_end: int
    dataitem: DataItem
class RFLP(DataValidationAbs):
    class DataItem(NamedTuple): seq_id: str; primer1_start: int; primer1_end: int; primer2_start: int; primer2_end: int; enzyme: str
    dataitem: DataItem
class PFGE(DataValidationAbs):
    class DataItem(NamedTuple): seq_id: str; restriction_enzyme: str
    dataitem: DataItem
class NOTE(DataValidationAbs):
    class DataItem(NamedTuple): note: str
    dataitem: DataItem


# fmt: on
class RecordEnum(Enum):
    # fmt: off
    SNP = SNP; SUB = SUB; DEL = DEL; INS = INS; MOB = MOB; AMP = AMP; CON = CON; INV = INV
    RA = RA; MC = MC; JC = JC; CN = CN; UN = UN
    CURA = CURA; FPOS = FPOS; PHYL = PHYL; TSEQ = TSEQ; PFLP = PFLP; RFLP = RFLP; PFGE = PFGE; NOTE = NOTE
    # fmt: on

    @classmethod
    def parse(
        cls,
        record_type: str,
        evidence_id: str | int,
        parent_ids: str | list[int] | None = ".",
        **kwargs,
    ):
        if isinstance(parent_ids, list):
            parent_ids = ",".join(str(i) for i in parent_ids)
        RE = cls[record_type].value
        pre = "\t".join((f"{kwargs[k]}" for k in RE.DataItem._fields))
        extra = "\t".join(
            f"{k}={v}" for k, v in kwargs.items() if k not in RE.DataItem._fields
        )
        return cls[record_type].value(evidence_id, parent_ids or "", f"{pre}\t{extra}")


Record = RecordEnum.parse


DATA2RECORD: dict[str, list[str]] = {}
for rtype, dtype in DataMeta.types.items():
    DATA2RECORD.setdefault(dtype, []).append(rtype)
del rtype, dtype


class RecordCollection:
    def __init__(self):
        self.SNP: list[SNP] = []
        self.SUB: list[SUB] = []
        self.DEL: list[DEL] = []
        self.INS: list[INS] = []
        self.MOB: list[MOB] = []
        self.AMP: list[AMP] = []
        self.CON: list[CON] = []
        self.INV: list[INV] = []
        self.RA: list[RA] = []
        self.MC: list[MC] = []
        self.JC: list[JC] = []
        self.CN: list[CN] = []
        self.UN: list[UN] = []
        self.CURA: list[CURA] = []
        self.FPOS: list[FPOS] = []
        self.PHYL: list[PHYL] = []
        self.TSEQ: list[TSEQ] = []
        self.PFLP: list[PFLP] = []
        self.RFLP: list[RFLP] = []
        self.PFGE: list[PFGE] = []
        self.NOTE: list[NOTE] = []

        self.index: dict[str | int | float, DataAbstract] = {}
        self.unindex: dict[str, list[DataAbstract]] = {}

    @classmethod
    def new(cls):
        return cls()  # type: ignore

    def set(self, record: DataAbstract):
        assert record.type in DataMeta.types
        getattr(self, record.type).append(record)
        if record.id in self.index:
            self.unindex.setdefault(record.id, []).append(self.index[record.id])
        if record.id in self.unindex or _convert_value(record.id) is None:
            self.unindex[record.id].append(record)
        else:
            self.index[record.id] = record

    def parents_of(self, record: DataAbstract | str | float | str):
        if not isinstance(record, DataAbstract):
            record = self[record]
            assert isinstance(record, DataAbstract)
        valid_pids = [
            pid for pid in record.parent_ids if _convert_value(record.id) is not None
        ]
        return [
            *(self.index[pid] for pid in valid_pids if pid in self.index),
            *(i for pid in valid_pids for i in self.unindex.get(pid, ())),
        ]

    @property
    def mutation(self) -> Iterable[DataMutationAbs]:
        return (i for rtype in DATA2RECORD["mutation"] for i in getattr(self, rtype))

    @property
    def evidence(self) -> Iterable[DataEvidenceAbs]:
        return (i for rtype in DATA2RECORD["evidence"] for i in getattr(self, rtype))

    @property
    def validation(self) -> Iterable[DataValidationAbs]:
        return (i for rtype in DATA2RECORD["validation"] for i in getattr(self, rtype))

    def __getitem__(self, item):
        try:
            return self.index[float(item)]
        except Exception:
            pass
        return self.index[item]

    def __len__(self):
        return sum(len(getattr(self, rtype)) for rtype in DataMeta.types)

    def __iter__(self):
        return (i for j in (self.mutation, self.evidence, self.validation) for i in j)
