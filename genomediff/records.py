import re

#########
## Fields should be kept synced with genomediff.cpp in the breseq source code

TYPE_SPECIFIC_FIELDS = {
    'SNP': ('seq_id', 'position', 'new_seq'),
    'SUB': ('seq_id', 'position', 'size', 'new_seq'),
    'DEL': ('seq_id', 'position', 'size'),
    'INS': ('seq_id', 'position', 'new_seq'),
    'MOB': ('seq_id', 'position', 'repeat_name', 'strand', 'duplication_size'),
    'AMP': ('seq_id', 'position', 'size', 'new_copy_number'),
    'CON': ('seq_id', 'position', 'size', 'region'),
    'INV': ('seq_id', 'position', 'size'),
    'RA': ('seq_id', 'position', 'insert_position', 'ref_base', 'new_base'),
    'MC': ('seq_id', 'start', 'end', 'start_range', 'end_range'),
    'JC': ('side_1_seq_id',
           'side_1_position',
           'side_1_strand',
           'side_2_seq_id',
           'side_2_position',
           'side_2_strand',
           'overlap'),
    'CN': ('seq_id', 'start', 'end', 'copy_number'),
    'UN': ('seq_id', 'start', 'end'),
    'CURA': ('expert',),
    'FPOS': ('expert',),
    'PHYL': ('gd',),
    'TSEQ': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end'),
    'PFLP': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end'),
    'RFLP': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end', 'enzyme'),
    'PFGE': ('seq_id', 'restriction_enzyme'),
    'NOTE': ('note',),
}

#########

class Metadata(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __repr__(self):
        return "Metadata({}, {})".format(repr(self.name), repr(self.value))

    def __str__(self, other):
        return "#={}\t{}".format(repr(self.name), repr(self.value))

    def __eq__(self, other):
        return self.__dict__ == other.__dict__



class Record(object):
    def __init__(self, type, id, document=None, parent_ids=None, **attributes):
        self.document = document
        self.type = type
        self.id = id
        self.parent_ids = parent_ids
        self.attributes = attributes

    @property
    def parents(self):
        if not self.parent_ids is None:
            return [self.document[pid] for pid in self.parent_ids]
        else:
            return []

    def __getattr__(self, item):
        try:
            return self.attributes[item]
        except KeyError:
            raise AttributeError


    def __repr__(self):
        return "Record('{}', {}, {}, {})".format(self.type,
                                             self.id,
                                             self.parent_ids,
                                             ', '.join('{}={}'.format(k, repr(v)) for k, v in self.attributes.items()))

    def __str__(self):

        parent_id_str='.'
        if not self.parent_ids is None:
          parent_id_str = ','.join(str(v) for v in self.parent_ids)

        ## Extract fields that are required (don't have field=value)
        remaining_items = self.attributes.copy()
        type_fields_key_list = TYPE_SPECIFIC_FIELDS[self.type]
        type_fields_list = []
        for type_field_key in type_fields_key_list:
          type_fields_list.append(str(self.attributes[type_field_key]))
          del remaining_items[type_field_key]
        type_fields_str = '\t'.join(type_fields_list)

        ## Remaining are printed as field=value
        attribute_str = '\t'.join('{}={}'.format(k, str(v)) for k, v in remaining_items.items())

        return '\t'.join([self.type,
                         str(self.id),
                         parent_id_str,
                         type_fields_str,
                         attribute_str
                         ])

    def __eq__(self, other):
        ''' this definition allows identical mutations in different genome diffs
            to be equal.'''
        return self.type == other.type and self.attributes == other.attributes

    def __ne__(self, other):
        return not self.__eq__(other)

    def satisfies(self, *args):
        '''
        Input: a variable number of conditions, e.g. 'gene_name==rrlA','frequency>=0.9'.
        Output: return true if all conditions are true (i.e. correspond to key-values in attributes.

        Find a condition that evaluates to false, otherwise return True.
        '''

        ## helper function to check if values are numbers
        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False

        for c in args:
            assert type(c) == str, "error: supplied condition is not a string."
            condition_pattern = re.compile(r'^(?P<key>[_a-z]+)'
                                            '(?P<comp>==|!=|<|<=|>|>=)'
                                            '(?P<val>[-_a-zA-Z0-9\.]+)')
            condition_match = condition_pattern.match(c)
            assert condition_match, "the supplied condition\n"+c+"\n could not be parsed."
            cond_key = condition_match.group('key')
            cond_comp = condition_match.group('comp')
            cond_val = condition_match.group('val')

            try: ## in case the given condition is not in the attributes.
                attribute_val = self.attributes[cond_key]
            except:
                continue

            ## add quote marks around strings before eval. can leave numbers alone.
            if not is_number(cond_val):
                cond_val = "\'"+cond_val+"\'"

            if not is_number(attribute_val):
                attribute_val = "\'"+attribute_val+"\'"
            else: ## attribute_val is a number in this record-- convert to str for eval.
                attribute_val = str(attribute_val)
            expr = attribute_val+cond_comp+cond_val
            if not eval(expr):
                return False
        return True
