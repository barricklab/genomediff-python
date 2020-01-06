from collections import OrderedDict
import re
from genomediff.records import TYPE_SPECIFIC_FIELDS, Metadata, Record

class GenomeDiffParser(object):
    def __init__(self, fsock=None, document=None):
        self._document = document
        self._fsock = fsock

    @staticmethod
    def _convert_value(value):
        for type_ in (int, float):
            try:
                return type_(value)
            except ValueError:
                pass

        if value == '.' or value == '':
            value = None
        return value

    def __iter__(self):
        metadata_pattern = re.compile(r'^#=(\w+)\s+(.*)$')
        mutation_pattern = re.compile(r'^(?P<type>[A-Z]{2,4})'
                                      '\t(?P<id>\d+)'
                                      '\t((?P<parent_ids>\d+(,\s*\d+)*)|\.?)'
                                      '\t(?P<extra>.+)?$')

        for i, line in enumerate(self._fsock):
            if not line:
                continue
            elif line.startswith('#'):
                match = metadata_pattern.match(line)
                if match:
                    yield Metadata(*match.group(1, 2))
            else:
                match = mutation_pattern.match(line)

                if match:
                    type = match.group('type')
                    id = int(match.group('id'))

                    parent_ids = match.group('parent_ids')
                    if parent_ids:
                        parent_ids = [int(id) for id in parent_ids.split(',')]

                    extra = match.group('extra').split('\t')
                    extra_dct = OrderedDict()

                    for name in TYPE_SPECIFIC_FIELDS[type]:
                        value = extra.pop(0)
                        extra_dct[name] = self._convert_value(value)

                    for k, v in (e.split('=', 1) for e in extra):
                        extra_dct[k] = self._convert_value(v)

                    yield Record(type, id, self._document, parent_ids, **extra_dct)
                else:
                    raise Exception('Could not parse line #{}: {}'.format(i, line))
