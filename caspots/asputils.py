import re

re_clause = re.compile("(\w+)\(([^\)]*)\)\.")
re_answer = re.compile("(\w+)\(([^\)]*)\)")

def parse_args(args):
    def parse_arg(a):
        a = a.strip()
        if a[0] == '"':
            return a.strip('"')
        else:
            return int(a)
    return [parse_arg(a) for a in args.split(",")]

