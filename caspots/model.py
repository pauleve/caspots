
from asputils import re_answer, parse_args

class Clause(set):
    def __str__(self):
        buf = " & ".join(["%sn_%s" % \
                ("!" if sign == -1 else "", var) for (var, sign) in self])
        if len(self) > 1:
            return "(%s)" % buf
        else:
            return buf

class DNF(list):
    def __init__(self, var):
        list.__init__(self)
        self.var = var
    def nusmv_expr(self):
        if len(self) == 0:
            return "FALSE"
        return " | ".join(map(str, self))
    def __str__(self):
        return "F_%s := %s;" % (self.var, self.nusmv_expr())

class BNModel:
    def __init__(self, id):
        self.id = id
        self.nodes = set()
        self.formula = {}
        self.clauses = {}
    def feed_from_asp(self, content):
        def key(e):
            try:
                return ["formula", "dnf", "clause"].index(e[0])
            except ValueError:
                return 3
        for (p, args) in sorted(content, key=key):
            args = parse_args(args)
            if p == "formula":
                node, fid = args
                self.nodes.add(node)
                self.formula[fid] = DNF(node)
            elif p == "dnf":
                fid, cid = args
                if cid not in self.clauses:
                    self.clauses[cid] = Clause()
                self.formula[fid].append(self.clauses[cid])
            elif p == "clause":
                cid, node, sign = args
                self.clauses[cid].add((node, sign))
            """
            elif p == "variable":
                self.nodes.add(args[0])
            """

    def default_false(self, nodes):
        fid = max(self.formula.keys())+1
        for node in nodes.difference(self.nodes):
            self.nodes.add(node)
            self.formula[fid] = DNF(node)
            fid += 1

    def __str__(self):
        buf = "%s %s %s\n" % ("#"*5, self.id, "#"*5)
        buf += "\n".join(map(str,self.formula.values()))
        return buf

def iter_bnmodels(inp, yield_line=False):
    bn = None
    for line in inp.readlines():
        if line.startswith("Answer:"):
            aid = int(line[7:])
            bn = BNModel(aid)
        elif bn is not None:
            bn.feed_from_asp(re_answer.findall(line))
            if yield_line:
                yield (bn, line)
            else:
                yield bn
            bn = None

if __name__ == "__main__":
    import sys
    for bn in iter_bnmodels(sys.stdin):
        print(bn)

