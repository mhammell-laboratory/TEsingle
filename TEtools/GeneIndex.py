import re
import gzip
import operator
import collections
import sys


class GFF_Reader():

    def __init__(self, filename, id_attribute):
        self.line_no = None
        self.filename = filename
        self.id_attribute = id_attribute
        self._re_attr_main = re.compile(r"\s*([^\s\=]+)[\s=]+(.*)")
        # this is a regular expression for chopping the description pairs into key and value

    def __iter__(self):
        self.line_no = 0
        if self.filename.lower().endswith((".gz", ".gzip")):
            lines = gzip.open(self.filename)
        else:
            lines = open(self.filename)

        for line in lines:
            self.line_no += 1
            if line == "\n" or line.startswith('#'):
                continue
            (seqname, source, feature, start, end, score, strand, frame, attributeStr) = line.split("\t")
            id = self.__parse_GFF_attr_string(attributeStr, self.id_attribute)

            yield id, seqname, strand, int(start), int(end), feature

        lines.close()
        self.line_no = None

    def __parse_GFF_attr_string(self, attributeStr, id_interested):
        for pairs in attributeStr.split(';'):
            if pairs.count('"') not in [0, 2]:
                raise ValueError("The attribute string seems to contain mismatched quotes.")
            nv = self._re_attr_main.match(pairs)  # YYY we could parse these differently and remove re dependence
            if not nv:
                raise ValueError("Failure parsing GFF attribute line.")
            val = nv.group(2)
            name = nv.group(1)
            if name == id_interested:
                return val
        return None

    def get_line_number_string(self):
        if self.line_no is None:
            return "file %s closed" % self.filename

        else:
            return "line %d of file %s" % (self.line_no, self.filename)


class IntervalTree(object):
    __slots__ = ('intervals', 'left', 'right', 'center')

    def __init__(self, intervals, depth=16, minbucket=16, _extent=None, maxbucket=512):

        depth -= 1
        if (depth == 0 or len(intervals) < minbucket) and len(intervals) < maxbucket:
            self.intervals = intervals
            self.left = self.right = None
            return

        if _extent is None:
            intervals.sort(key=operator.attrgetter('start'))

        left, right = _extent or (intervals[0].start, max(i.stop for i in intervals))
        center = (left + right) / 2.0

        self.intervals = []
        lefts, rights = [], []

        for interval in intervals:
            if interval.stop < center:
                lefts.append(interval)
            elif interval.start > center:
                rights.append(interval)
            else:  # overlapping.
                self.intervals.append(interval)

        self.left = lefts and IntervalTree(lefts, depth, minbucket, (intervals[0].start, center)) or None
        self.right = rights and IntervalTree(rights, depth, minbucket, (center, right)) or None
        self.center = center

    def find(self, start, stop):
        if self.intervals and not stop < self.intervals[0].start:
            overlapping = [i for i in self.intervals if i.stop >= start and i.start <= stop]
        else:
            overlapping = []

        if self.left and start <= self.center:
            overlapping += self.left.find(start, stop)

        if self.right and stop >= self.center:
            overlapping += self.right.find(start, stop)

        return overlapping

    def find_gene(self, start, stop):
        if self.intervals and not stop < self.intervals[0].start:
            overlapping = [i.gene for i in self.intervals if i.stop >= start and i.start <= stop]
        else:
            overlapping = []

        if self.left and start <= self.center:
            overlapping += self.left.find_gene(start, stop)

        if self.right and stop >= self.center:
            overlapping += self.right.find_gene(start, stop)

        return overlapping

    def __iter__(self):
        if self.left:
            for l in self.left:
                yield l

        for i in self.intervals:
            yield i

        if self.right:
            for r in self.right:
                yield r


class Interval(object):
    def __init__(self, gene_id, start, stop):
        self.start = start
        self.stop = stop
        self.gene = gene_id

    def __repr__(self):
        return "%s (%i, %i)" % (self.gene, self.start, self.stop)


class GeneFeatures:
    def __init__(self, GTFfilename, stranded, feature_type, id_attribute):

        self.featureIdxs_plus = {}
        self.featureIdxs_minus = {}
        self.featureIdxs_nostrand = {}
        self.features = []
        self.lengths = []
        self.featureDict = {}
        self.lengthDict = {}

        self.read_features(GTFfilename, stranded, feature_type, id_attribute)

    def read_features(self, gff_filename, stranded, feature_type, id_attribute):

        temp_plus = collections.defaultdict(dict)
        temp_minus = collections.defaultdict(dict)
        temp_nostrand = collections.defaultdict(dict)

        gff = GFF_Reader(gff_filename, id_attribute)
        i = 0
        counts = 0
        try:
            for f in gff:
                if f[0] is None:
                    continue
                if f[5] == feature_type:
                    counts += 1
                    if stranded != "no" and f[2] == ".":
                        sys.stderr.write("Feature %s does not have strand information." % (f[0]))
                    try:
                        if f[2] == ".":
                            temp_nostrand[f[1]][f[0]].append((f[3], f[4]))
                    except:
                        temp_nostrand[f[1]][f[0]] = [(f[3], f[4])]

                    try:
                        if f[2] == "+":
                            temp_plus[f[1]][f[0]].append((f[3], f[4]))
                    except:
                        temp_plus[f[1]][f[0]] = [(f[3], f[4])]

                    try:
                        if f[2] == "-":
                            temp_minus[f[1]][f[0]].append((f[3], f[4]))
                    except KeyError:
                        temp_minus[f[1]][f[0]] = [(f[3], f[4])]

                    # save gene id
                    if f[0] not in self.features:
                        self.features.append(f[0])
                    if f[0] not in self.lengthDict:
                        self.lengthDict[f[0]] = 0
                    self.lengthDict[f[0]] += f[4] - f[3]

                    i += 1
                    if i % 100000 == 0:
                        print("%d GTF lines processed.\n" % i, flush=True)
        except:
            sys.stderr.write("Error occured in %s.\n" % gff.get_line_number_string())
            raise

        if counts == 0:
            sys.stderr.write("Warning: No features of type '%s' found in gene GTF file.\n" % feature_type)

        # build interval trees
        for each_chrom in temp_plus:
            inputlist = []
            for each_gene in temp_plus[each_chrom]:
                for (start, end) in temp_plus[each_chrom][each_gene]:
                    inputlist.append(Interval(each_gene, start, end))
            self.featureIdxs_plus[each_chrom] = IntervalTree(inputlist)

        for each_chrom in temp_minus:
            inputlist = []
            for each_gene in temp_minus[each_chrom]:
                for (start, end) in temp_minus[each_chrom][each_gene]:
                    inputlist.append(Interval(each_gene, start, end))
            self.featureIdxs_minus[each_chrom] = IntervalTree(inputlist)

        for each_chrom in temp_nostrand:
            inputlist = []
            for each_gene in temp_nostrand[each_chrom]:
                for (start, end) in temp_nostrand[each_chrom][each_gene]:
                    inputlist.append(Interval(each_gene, start, end))
            self.featureIdxs_nostrand[each_chrom] = IntervalTree(inputlist)

        self.featureDict = dict(zip(self.features, range(len(self.features))))
        for elem in self.featureDict:
            self.lengths.append(self.lengthDict[elem])

    def getFeatures(self):
        return self.features

    def getLengths(self):
        return self.lengths

    def Gene_annotation(self, itv_list):
        genes = []
        for itv in itv_list:
            fs = []
            try:
                if itv[3] == "+":
                    if itv[0] in self.featureIdxs_plus:
                        fs = self.featureIdxs_plus[itv[0]].find_gene(itv[1], itv[2])

                if itv[3] == "-":
                    if itv[0] in self.featureIdxs_minus:
                        fs = self.featureIdxs_minus[itv[0]].find_gene(itv[1], itv[2])

                if itv[3] == ".":
                    if itv[0] in self.featureIdxs_minus:
                        fs = self.featureIdxs_minus[itv[0]].find_gene(itv[1], itv[2])

                    if itv[0] in self.featureIdxs_plus:
                        fs += self.featureIdxs_plus[itv[0]].find_gene(itv[1], itv[2])
                    if itv[0] in self.featureIdxs_nostrand:
                        fs += self.featureIdxs_nostrand[itv[0]].find_gene(itv[1], itv[2])

                if len(fs) > 0:
                    genes = genes + fs

            except:
                raise

        gen_vect = []
        for gen in set(genes):
            gen_vect.append(self.featureDict[gen])

        return gen_vect
