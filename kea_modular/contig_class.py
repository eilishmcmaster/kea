
def contig_class_maker(ref_contig, nucmer_array):

    class contig(ref_contig, nucmer_array):
        def __init__(self,name):
            self.name = name

        def reference(self, name, ref_contig):
            if name is ref_contig:
                self.reference = True
            else:
                self.reference = False


        start = 1
        def end(self, ref_contig):

        matches = []
        def add_match(self, match):
            self.matches.append(match)

