from pyparsing import Word, OneOrMore, Group, ParseException


##### pyparsing based cigar parser
def convertIntegers(tokens):
	return int(tokens[0])


alt_type_tokens = "SIMDX"
digits = "0123456789"
an_alt = Word(alt_type_tokens)
alt_length = Word(digits).setParseAction(convertIntegers)
alt_and_length = Group(alt_length + an_alt)
cigar_string_parser = OneOrMore(alt_and_length)


def parse_cigar_string(cigar_string):
	return cigar_string_parser.parseString(cigar_string)



# Exemple usage:
# cigar_string.parseString("76M1I257M1I22M1I99M")
# % timeit cigar_string.parseString("76M1I257M1I22M1I99M")[1]
