#!/usr/bin/env python3

"""
generate_mappings.py
reads over BLAST results then gets the top four hits per read, accounting for
direction and reading frames (i.e., up to 2 hits per direction per read)
"""

class Blast():
    def __init__(self, line):
        self.line = line.strip()
        
        fields = self.line.split()
        self.query=fields[0]
        self.subject=fields[1]
        self.sstart=int(fields[9])
        self.send=int(fields[10])
        self.bitscore=float(fields[12])
        self.evalue=float(fields[11])
        self.line=line

    def overlap(self, other):
        """returns the amount that the results overlap by on the subject"""
        range1 = (self.sstart, self.ssend)
        range2 = (other.sstart, other.ssend)
        return max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]))

    def is_better(self, other, by="evalue"):
        """true if the evalue/score of self is better than other"""
        if by == evalue:
            return self.evalue < other.evalue
        elif by == bitscore:
            return self.bitscore > other.bitscore
        else:
            raise ValueError("\'by\' must be either \'evalue\' or \'bitscore\'")

    def frame(self):
        """returns the frame number of the blast query"""
        return int(self.query[-1])


class Read:
    def __init__(self):
        self.forward = []   # List with at most two elements
        self.reverse = []   # List with at most two elements

    def _insert_forward(self, hit):
        counter = 0
        if len(self.forward == 2):
            if hit.is_better(self.forward[0]):
                self.forward[1] = self.forward[0]
                self.forward[0] = hit
            elif hit.is_better(self.forward[1]):
                self.forward[1] = hit
        elif len(self.foward == 1):
            if hit.is_better(self.forward[0]):
                temp = self.forward[0]
                self.forward[0] = hit
                self.forward.append(temp)
        elif len(self.forward == 0):
            self.forward.append(hit)

    def _insert_reverse(self, hit):
        counter = 0
        if len(self.reverse == 2):
            if hit.is_better(self.reverse[0]):
                self.reverse[1] = self.reverse[0]
                self.reverse[0] = hit
            elif hit.is_better(self.reverse[1]):
                self.reverse[1] = hit
        elif len(self.reverse == 1):
            if hit.is_better(self.reverse[0]):
                temp = self.reverse[0]
                self.reverse[0] = hit
                self.reverse.append(temp)
        elif len(self.reverse == 0):
            self.reverse.append(hit)

    def insert(self, blast_line):
        hit = Blast(blast_line)
        if hit.frame() < 4:
            self._insert_forward(hit)
        else:
            self._insert_reverse(reverse)      



