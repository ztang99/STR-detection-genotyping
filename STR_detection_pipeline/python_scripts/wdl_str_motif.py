class STRMotif:
    def __init__(self, gene=None, motif=None, start=None, end=None, chrom=None, carriers=None):
        self.gene = gene
        self.motif = motif 
        self.start = start
        self.end = end
        self.chrom = chrom
        self.carriers = carriers if carriers else []
        
    @classmethod
    def from_bed_line(cls, line):
        """Create STRMotif from a bed file line"""
        chrom, start, end, gene, motif = line.strip().split('\t')
        return cls(
            gene=gene,
            motif=motif,
            start=int(start),
            end=int(end),
            chrom=chrom
        )
    
    @classmethod
    def from_ehdn_eh_results(cls, ehdn_eh_combined):
        """Create STRMotif from combined EHdn and EH results"""
        return cls(
            gene=ehdn_row['gene'],
            motif=ehdn_row['motif'],
            start=ehdn_row['start'],
            end=ehdn_row['end'],
            chrom=ehdn_row['chr']
        )
    
    def add_carrier(self, bam_path):
        """Add a carrier (bam file path) to this motif"""
        if bam_path not in self.carriers:
            self.carriers.append(bam_path)
    
    def to_dict(self):
        """Convert to dictionary for easier serialization"""
        return {
            'gene': self.gene,
            'motif': self.motif,
            'start': self.start,
            'end': self.end,
            'chrom': self.chrom,
            'carriers': self.carriers
        }
    
    @classmethod
    def from_dict(cls, d):
        """Create STRMotif from dictionary"""
        return cls(**d)
    
    def get_region(self):
        """Get region in format chr:start-end"""
        return f"{self.chrom}:{self.start}-{self.end}"
    
    def __str__(self):
        return f"{self.gene}_{self.motif}_{self.chrom}:{self.start}-{self.end}"