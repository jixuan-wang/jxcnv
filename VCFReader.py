import vcf

class VCFReader:
    def __init__(self, vcf_file):
        try:
            vfile = open(vcf_file, 'r')
            self.vcf_reader = vcf.Reader(vfile)
        except IOException:
            print 'This VCF file does not exist: ' + vcf_file

    def chrInt2Str(self, chromosome_int):
        try:
            if int(chromosome_int) == 23:
                return 'X'
            elif int(chromosome_int) == 24:
                return 'Y' 
            else:
                return str(chromosome_int)
        except:
            return str(chromosome_int)
        

    # Given a region, check the number of sites which are heterozygous calls
    def getSNPInfo(self, sample_id, intervals_list):
        snp_info = list()
        if self.vcf_reader:
            for in_list in intervals_list:
                for interval in in_list:
                    # The start and end coordinates are in the zero-based, half-open coordinate system
                    records = self.vcf_reader.fetch(self.chrInt2Str(interval._chr), interval._bp1, interval._bp2 + 1) 
                    if records:
                        het_num = 0
                        for record in records:
                            genotype = record.genotype(sample_id)
                            if genotype and genotype.is_het:
                                het_num += 1    
                        snp_info.append(het_num)
                    else:
                        snp_info.append(0)
                        
        return snp_info
