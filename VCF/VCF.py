import os

from typing import Dict, Generator, List, Literal, Type

from .varRecord import varRecord
import gzip
import pandas as pd


__all__ = ('metaFILTER', 'metaFORMAT', 'metaINFO', 'VCF')

class metaFILTER:
    """Class contains 'FILTER' field meta information
    
    Example
    -------
    ##FILTER=<ID=PASS,Description="All filters passed">
    """
    
    def __init__(self,
        id : str,
        desc : str):
        
        self.id = id
        self.desc = desc
        
        return
    
    def __str__(self) -> str:
        return f'<ID={self.id},Description={self.desc}>'


class metaFORMAT:
    """Class contains 'FORMAT' field meta information
    
    Possible Types of FORMAT : Integer, Float, Character, String
    
    Example
    -------
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic ..">
    """
    
    TYPES = {'Integer', 'Float', 'Character', 'String'}
    
    def __init__(self,
        id : str,
        number : str,
        type : str,
        desc : str):

        self.id = id
        self.number = number
        self.type = type
        self.desc = desc

        return
    
    @property
    def type(self) -> str:
        return self._type
    
    @type.setter
    def type(self, _type: str):
        if _type not in self.TYPES:
            print(f'[Warning] FORMAT field should have 4 types ({self.TYPES})')
            print(f'This line contains {_type}')
        self._type = _type
        return

    def __str__(self) -> str:
        return f'<ID={self.id},Number={self.number},Type={self.type},Description={self.desc}>'


class metaINFO:
    """Class contains 'INFO' field meta information
    
    Possible Types of INFO : Integer, Float, Flag, Character, String
    
    Example
    -------
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate...">
    """

    TYPES = {'Integer', 'Float', 'Flag', 'Character', 'String'}


    def __init__(self,
        id : str,
        number : str,
        type : str,
        desc : str):
        
        self.id = id
        self.number = number
        self.type = type
        self.desc = desc
        
        return
    
    @property
    def type(self) -> str:
        return self._type
    
    @type.setter
    def type(self, _type: str):
        if _type not in self.TYPES:
            print(f'[Warning] INFO field should have 5 types ({self.TYPES})')
            print(f'This line contains {_type}')
        self._type = _type
        return
    
    def __str__(self) -> str:
        return f'<ID={self.id},Number={self.number},Type={self.type},Description={self.desc}>'
    

class VCF:
    """This supports various functions for handling VCF file.
    VCF file contains variant information from sample (or samples),
    each line represent each genomic variant.
    VCF consist of meta information lines, header line, and data lines.
    Meta information lines start with '##' and header line start with '#'.
    (Detail explanation of meta lines are below)
    Data lines which follow header line do not start with '#', and 
    each data line actually represent each variant.
    Header line and data lines are tab seperated and each column
    is called as field. The name of columns(fields) are determined at
    header line. For example,
    
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO    FORMAT  SAMPLE1...
    
    The first 8 fields are mandatory.
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
    
    And most of VCFs produced by normal analysis pipeline have fields
    'FORMAT', 'SAMLE1', 'SAMLE2', ...
    
    Of 9 fields, 'FILTER', 'INFO', and 'FORMAT' fields contain
    more complicated information. For example,
    
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    Y	2728456	rs2058276	T	C	32	.	AC=2;AN=2;DB;DP=182;H2;NS=65
    Y	2734240	.	G	A	31	.	AC=1;AN=2;DP=196;NS=63
    Y	2743242	.	C	T	25	.	AC=1;AN=2;DP=275;NS=66
    Y	2746727	.	A	G	34	.	AC=2;AN=2;DP=179;NS=64
    Y	2777970	.	T	A	67	.	AC=1;AN=2;DP=225;NS=67

    'INFO' field contain a few of '[key]=[value]' pairs where
    each key indicates statistics of a variant. Details for each
    key are represented at the top of VCF, lines starting with '##'.
    These lines contain meta information for data line of VCF, and
    are written at the top of file. For example,
    
    ##fileformat=VCFv4.0
    ##fileDate=20100610 
    ##source=glfTools v3
    ##reference=1000GenomesPilot-NCBI36 
    ##phasing=NA
    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Mapped Reads">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
    ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
    ##FILTER=<ID=NUYR,Description="Variant in non-unique Y region">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype	Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
    ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    
    The format of meta lines should be [key]=[value], and
    [value] could also have [key]=[value] pairs in '<>'.
    
    Meta lines are important as data lines, BI.VCF.VCF class supports
    method which parses meta lines and returns dictionary.
    
    Example
    -------
    >>> from BI import VCF
    >>> vcf = VCF.VCF('./data/small.vcf')
    >>> vcf.sanity_check()
    AC of INFO has unknown range
    True
    >>> vcf.open()
    >>> meta_info = vcf.parse_meta_info_lines()
    >>> for key in sorted(list(meta_info.keys())):
    ...     if key in {'FILTER', 'INFO', 'FORMAT'}:
    ...         for id in sorted(list(meta_info[key].keys())):
    ...             print(key, id, meta_info[key][id].desc)
    FILTER NUYR "Variant in non-unique Y region"
    FORMAT DP "Depth"
    FORMAT GQ "Genotype	Quality"
    FORMAT GT "Genotype"
    INFO AC "Allele count in genotypes"
    INFO AN "Total number of alleles in called genotypes"
    INFO DB "dbSNP membership, build 129"
    INFO DP "Total Depth"
    INFO H2 "HapMap2 membership"
    INFO NS "Number of Samples With Mapped Reads"
    >>> for line in vcf.reader():
    ...     print(line)
    varRecord(Y, 2728456, rs2058276, T, C, 32.0, ., {'AC': '2', 'AN': '2', 'DB': '', 'DP': '182', 'H2': '', 'NS': '65'})
    varRecord(Y, 2734240, ., G, A, 31.0, ., {'AC': '1', 'AN': '2', 'DP': '196', 'NS': '63'})
    varRecord(Y, 2743242, ., C, T, 25.0, ., {'AC': '1', 'AN': '2', 'DP': '275', 'NS': '66'})
    varRecord(Y, 2746727, ., A, G, 34.0, ., {'AC': '2', 'AN': '2', 'DP': '179', 'NS': '64'})
    varRecord(Y, 2777970, ., T, A, 67.0, ., {'AC': '1', 'AN': '2', 'DP': '225', 'NS': '67'})
    varRecord(Y, 2782506, rs2075640, A, G, 38.0, ., {'AC': '1', 'AN': '2', 'DB': '', 'DP': '254', 'H2': '', 'NS': '66'})
    varRecord(Y, 2783755, ., G, A, 51.0, ., {'AC': '1', 'AN': '2', 'DP': '217', 'NS': '67'})
    varRecord(Y, 2788927, rs56004558, A, G, 38.0, ., {'AC': '1', 'AN': '2', 'DB': '', 'DP': '173', 'NS': '60'})
    >>> df = vcf.load_to_dataframe()
    """

    def __init__(self, vcf_: str) -> None:
        """Initialize VCFProcessor class
        It sets default self.header, checks compressness, opens vcf.
        It prepares vcf processing jobs.

        Parameters
        ----------
        vcf_ : str
            Absolute path of vcf for processing
        """

        ## default header
        self.header = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                    'QUAL', 'FILTER', 'INFO', 'FORMAT',
                    'SAMPLE']
        self.vcf = vcf_
        self._compressed = self._chk_compressed(self.vcf)
        self.f_obj = False
        self.mode = False

        return

    @property
    def vcf(self) -> str:
        return self._vcf

    @vcf.setter
    def vcf(self, path: str) -> None:
        ## check existence
        if os.path.isfile(path):
            self._vcf = path
        else:
            print(f'[WARNING] {path} does not exist...')
            self._vcf = path
        return

    def _chk_compressed(self, vcf_: str) -> bool:
        """Check whether vcf is compressed (by extension)

        Parameters
        ----------
        vcf_ : str
            Path of input vcf

        Returns
        -------
        bool
            True if compressed, else False
        """

        extension = vcf_.split('.')[-1]
        if extension == 'vcf':
            return False
        elif extension == 'gz':
            return True
        else:
            raise Exception('Check file extension (only .vcf | .vcf.gz)')

    def sanity_check(self) -> bool:
        """Check integrity of self.vcf file.
        Checking codes are from the content of official VCF format docs
        (https://samtools.github.io/hts-specs/VCFv4.2.pdf).

        > Load ## header lines and check it's ID and Number... etc. is normally matched
         > Fields of INFO/FORMAT are in meta info lines (V, only INFO)
         > In INFO/FORMAT, Type must be Integer, Float, Flag, Character, String
         > In INFO/FORMAT, Number must be 0, 1, ..., A, R, G, .
         > If Type is Flag -> Number must be 0
        > All variant lines have the same # of columns (V)
        > Whether REF == ALT allele (V)
        > Check DP == 0 (V)
        
        Returns
        -------
        bool
            Is self.vcf normal format
        """
        
        meta_info = self.parse_meta_info_lines()
        if not meta_info:
            print(f'[Warning] VCF does not have meta info lines')
            return False
        info_fields = set([id for id in meta_info['INFO'].keys()])
        
        ## check meta information lines
        for field in meta_info['INFO'].values():
            ## check Number == '.' (unknown value)
            if field.number == '.':
                print(f'{field.id} of INFO has unknown range')
            ## check strange Type
            if field.type not in metaINFO.TYPES:
                print(f'{field.id} of INFO is unknown type')
                print(f'{field.type}')
        if 'FORMAT' in meta_info.keys():
            for field in meta_info['FORMAT'].values():
                ## check Number == '.'
                if field.number == '.':
                    print(f'{field.id} of FORMAT has unknown range')
                ## check strange Type
                if field.type not in metaFORMAT.TYPES:
                    print(f'{field.id} of FORMAT is unknown type')
                    print(f'{field.type}')
        
        vcf = VCF(self.vcf)
        vcf.open()
        for variant in vcf.reader():
            ## check column count
            if len(vcf.header) != variant.get_column_num():
                print(f'[Warning] VCF has a line with different column count')
                print(f'{variant}')
                return False
            ## check whether REF == ALT
            if variant.ref == variant.alt:
                print(f'[Warning] VCF has a line with ref == alt')
                print(f'{variant}')
                return False
            ## check if DP == 0 (in INFO)
            if 'DP' in variant.info.keys():
                if variant.info['DP'] == 0:
                    print(f'[Warning] VCF has a line with DP == 0')
                    print(f'{variant}')
                    return False
            ## check INFO fields and meta lines
            for key in variant.info.keys():
                if key not in info_fields:
                    print(f'[Warning] Detect strange INFO field')
                    print(f'Check {key} is in meta lines')
                    print(f'{variant}')
                    return False
                
        return True
    
    # def parse_meta_info_lines(self) -> Dict[str, List]:
    def parse_meta_info_lines(self) -> Dict[str, dict] or Dict[str, list]:
        """Parse meta information lines starting with '##' to save meta info of VCF file.
        Meta information lines must be key=value pairs. Please refer to official VCF format
        document.
        
        Meta information lines example,
        
        ##fileformat=VCFv4.0
        ##fileDate=20100610 
        ##source=glfTools v3
        ##reference=1000GenomesPilot-NCBI36 
        ##phasing=NA
        ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Mapped Reads">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
        ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
        ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
        ##FILTER=<ID=NUYR,Description="Variant in non-unique Y region">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype	Quality">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
        ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        
        Returns
        -------
        Dict[str, dict]
            Dictionary containing meta information
        """
        
        meta_info = dict()
        
        ## init new VCF instance
        vcf = VCF(self.vcf)
        vcf.open()
        line = vcf.readline(skip_header=False)
        while line != '':
            if line[:2] != '##':
                break
            
            ## '##FILTER=<...>' -> 'FILTER'
            field = line[2:].strip().split('=')[0]
            if field not in meta_info.keys():
                # meta_info[field] = set()
                # meta_info[field] = []
                if field in {'FILTER', 'INFO', 'FORMAT'}:
                    meta_info[field] = dict()
                else:
                    meta_info[field] = []
            
            ## '##FILTER=<ID=..,Description="..">' -> '<ID=..,Description="..">'
            contents: str = '='.join(line[2:].strip().split('=')[1:])
            if contents[0] == '<':
                ## '<ID=..,Description=".."">' -> ['ID=..', 'Description=".."']
                temp = contents[1:-1].split(',')
                
                ## 'Description' can have ','
                ## Ex) Description="dbSNP membership, build 129"
                contents_list = []
                for content in temp:
                    if '=' in content:
                        contents_list.append(content)
                    else:
                        contents_list[-1] += f',{content}'

                temp = dict()
                for content in contents_list:
                    # print(content)
                    key = content.split('=')[0]
                    value = content.split('=')[1]
                    ## {'ID':'..', 'Description':'..'}
                    temp[key] = value
                
                ## '##FILTER=<ID=..,Description=..>'
                ## -> {'FILTER': [metaFILTER, ...]}
                if field == 'FILTER':
                    # meta_info[field].append(metaFILTER(temp['ID'], temp['Description']))
                    meta_info[field][temp['ID']] = metaFILTER(temp['ID'], temp['Description'])
                elif field == 'FORMAT':
                    # meta_info[field].append(metaFORMAT(temp['ID'], temp['Number'],\
                    #                         temp['Type'], temp['Description']))
                    meta_info[field][temp['ID']] = metaFORMAT(temp['ID'], temp['Number'],\
                                            temp['Type'], temp['Description'])
                elif field == 'INFO':
                    # meta_info[field].append(metaINFO(temp['ID'], temp['Number'],\
                    #                         temp['Type'], temp['Description']))
                    meta_info[field][temp['ID']] = metaINFO(temp['ID'], temp['Number'],\
                                            temp['Type'], temp['Description'])
                else:
                    ## '##GATKCommandLine=<ID=..,CommandLine=..>'
                    ## -> {'GATKCommandLine': [{'ID':'..', 'CommandLine':'..'}]}
                    meta_info[field].append(temp)
            else:
                ## '##VCFformat=4.2v' -> {'VCFformat':['4.2v']}
                meta_info[field].append(contents)
            line = vcf.readline(skip_header=False)
        
        return meta_info
    
    def reader(self) -> Generator[varRecord, None, None]:
        """Generator function read vcf file line by line.
        
        Yields
        ------
        varRecord
            varRecord parsed from each vcf line
        """
        
        if not self.f_obj:
            print(f'[ERROR] {self.vcf} is not opened.')
            return
        
        meta_info = self.parse_meta_info_lines()

        line = self.readline()
        while line != '':
            cols = line.strip().split('\t')
            if len(cols) < 9:
                yield varRecord(*cols[:8], meta_info)
            else:
                yield varRecord(*cols[:8], meta_info, cols[8], cols[9:], self.header[9:])
            line = self.readline()

        return

    def open(self,
        mode : Literal['r', 'w'] = 'r'
        ) -> None:
        
        """Open self.vcf (.vcf or .vcf.gz path) and
        save TextIOWrapper instance by attribute (self.f_obj).

        Parameters
        ----------
        mode: str (default: 'r')
            'r': read, 'w': write
        """
        
        if mode not in {'r', 'w'}:
            raise Exception('Only \'r\' and \'w\' mode supports')

        ## .vcf.gz case
        if self._compressed:
            mode = 'rb' if mode == 'r' else 'wb'
            f_obj = gzip.open(self.vcf, mode)
            self.mode = mode
        else:
            f_obj = open(self.vcf, mode)
            self.mode = mode
        self.f_obj = f_obj
        return

    def close(self) -> None:
        """Close self.vcf"""
        
        self.f_obj.close()
        self.f_obj = False
        
        return

    def readline(self, skip_header=True) -> str:
        """Read vcf file by one line (~= TextIOWrapper.readline())

        Parameters
        ----------
        skip_header : bool, optional
            Boolean val whether skip header line which starts with '#', by default True

        Returns
        -------
        str
            A line from vcf
        """

        if not self.f_obj:
            self.open(mode='r')

        if self._compressed:
            line = self.f_obj.readline().decode()
            if skip_header:
                while line[:2] == '##':
                    line = self.f_obj.readline().decode()
        else:
            line = self.f_obj.readline()
            if skip_header:
                while line[:2] == '##':
                    line = self.f_obj.readline()

        if line == '':
            return line
        
        if line[0] == '#' and line[:2] != '##':
            self.header = line[1:].strip('\n').split('\t')
            if skip_header:
                line = self.f_obj.readline().decode() if self._compressed else self.f_obj.readline()

        return line

    def write(self, line : str) -> None:
        """Write an input line at self.vcf file.
        This method check input line whether it contains
        proper columns for vcf format comparing to 'self.header'.
        
        Parameters
        ----------
        line : str
            A line for writing vcf file
        """

        ## check input line if it has proper format
        if line[:2] != '##':
            if line[0] == '#':
                self.header = line[1:].strip('\n').split('\t')
            else:
                cols = line.strip('\n').split('\t')
                if len(cols) != len(self.header):
                    print('Skip improper formatted line')
                    return
        
        if not self.f_obj:
            self.open(mode='w')

        if line[-1] != '\n':
            line += '\n'
        
        if self._compressed:
            self.f_obj.write(line.encode())
        else:
            self.f_obj.write(line)

        return

    def is_indel(self, line: str) -> bool:
        """Return whether a variant(line) is indel or not"""

        A_allele_idx = self.header.index('REF')
        B_allele_idx = self.header.index('ALT')

        cols     = line.strip('\n').split('\t')
        A_allele = cols[A_allele_idx]
        B_allele = cols[B_allele_idx]

        return (len(A_allele) != len(B_allele))

    def is_snp(self, line: str) -> bool:
        """Return whether a variant(line) is snp or not"""

        A_allele_idx = self.header.index('REF')
        B_allele_idx = self.header.index('ALT')

        cols     = line.strip('\n').split('\t')
        A_allele = cols[A_allele_idx]
        B_allele = cols[B_allele_idx]

        return (len(A_allele) == 1 and len(B_allele) == 1)

    def is_mnp(self, line: str) -> bool:
        """Return whether a variant(line) is mnp or not (snp included)"""

        A_allele_idx = self.header.index('REF')
        B_allele_idx = self.header.index('ALT')

        cols = line.strip('\n').split('\t')
        A_allele = cols[A_allele_idx]
        B_allele = cols[B_allele_idx]

        return (len(A_allele) == len(B_allele))
        
    def get_genotype(self,
        genotypes : List[str]
        ) -> Generator[str, None, None]:
        
        """Bring a variant line containing input genotype in 'GT' field.
        
        Parameters
        ----------
        genotype : List[str]
            Input genotype list for checking
        
        Yields
        ------
        Generator[str, None, None]
            Line containing input genotype        
        """
        
        if not self.f_obj:
            ## Open VCF first
            return
        
        line = self.readline()
        if 'FORMAT' not in self.header:
            ## GT info not in VCF
            return
        
        format_idx = self.header.index('FORMAT')
        while line != '':
            cols = line.strip('\n').split('\t')
            format = cols[format_idx].split(':')
            try:
                gt_idx = format.index('GT')
            except ValueError:
                ## GT info not in VCF
                return
            
            chk_gt = False
            for spl_idx in range(format_idx+1, len(self.header)):
                spl_genotype = cols[spl_idx].split(':')[gt_idx]
                if spl_genotype in genotypes:
                    chk_gt = True
            if chk_gt:
                yield line
            line = self.readline()
        
        return

    def _is_header(self, line : str) -> bool:
        """Check whether input line is vcf header line.
        
        Parameters
        ----------
        line : str
            Input vcf line for checking
        
        Returns
        -------
        bool
            Whether line is header line
        """
        
        return line[0] == '#'

    def _set_header(self) -> None:
        """Set header of VCF file"""
        
        vcf = VCF(self.vcf)
        vcf.open()
        line = vcf.readline(skip_header=False)
        while line != '':
            if line[:2] == '##':
                line = vcf.readline(skip_header=False)
                continue
            if line[0] == '#':
                break
            line = vcf.readline(skip_header=False)
        vcf.close()
        
        header = line[1:].strip().split('\t')
        self.header = header
        
        return
    
    def load_to_dataframe(self) -> Type['pd.DataFrame']:
        """Load variant information to pandas dataframe
        * consider memory usage for large vcf file
        
        Test
        ----
        >>> from BI import VCF
        >>> vcf = VCF.VCF('./data/small.vcf')
        >>> df = vcf.load_to_dataframe()
        >>> df
        """
        
        vcf = VCF(self.vcf)
        vcf.open()
        
        ## load meta information
        meta_info = vcf.parse_meta_info_lines()
        
        vcf._set_header()
        df_columns = []
        for col in vcf.header[:8]:
            if col == 'INFO':
                for info_id in sorted(list(meta_info[col].keys())):
                    ## INFO -> INFO-NS, INFO-DP,...
                    df_columns.append(f'{col}-{info_id}')
            else:
                df_columns.append(col)

        ## optional fields
        if len(vcf.header) > 8:
            format_columns = []
            for format_id in sorted(list(meta_info['FORMAT'].keys())):
                format_columns.append(f'{format_id}')
                
            for sample in vcf.header[9:]:
                for format in format_columns:
                    df_columns.append(f'{sample}-{format}')

        list_to_df = []
        for var in vcf.reader():
            ## mandatory fields
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            chrom = var.chrom
            pos   = var.pos
            id    = var.ID
            ref   = var.ref
            alt   = var.alt
            qual  = var.qual
            filt  = ';'.join(var.filter)
            
            ## parsing 'INFO' values
            info_vals = []
            info_cols = [ col for col in df_columns if col.startswith('INFO-') ]
            for col in info_cols:
                try:
                    info_val = var.info[f'{"-".join(col.split("-")[1:])}']
                except KeyError:
                    info_val = 'X'
                info_vals.append(info_val)
                
            variant = [chrom, pos, id, ref, alt, qual, filt, *info_vals]
            
            ## optional fields
            #FORMAT  SAMPLE1 ...
            if len(vcf.header) > 8:
                format_cols = [ col for col in df_columns if col.startswith('FORMAT-') ]
                for sample in vcf.header[9:]:
                    sample_cols = []
                    for col in format_cols:
                        try:
                            sample_val = var.sample_info[sample][f'{"-".join(col.split("-")[1:])}']
                        except KeyError:
                            sample_val = 'X'
                        sample_cols.append(sample_val)
                
                    variant += sample_cols
            list_to_df.append(variant)
        df = pd.DataFrame(data=list_to_df, columns=df_columns)
        
        return df


def _test():
    """test by CI script"""
    
    import doctest
    doctest.testmod()
    
    return


if __name__ == '__main__':
    _test()