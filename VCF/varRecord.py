"""
This is for loading tab-seperated variant information
in VCF file. The specification of VCF format follows VCF v4.3.
Each line in VCF has multiple fields, first 8 columns are
mandatory fields and other columns are optional fields.
"""

import sys
from typing import List, Dict, TypeVar

__all__ = ('varRecord',)

INFO = TypeVar('INFO', int, float, str)
FORMAT = TypeVar('FORMAT', int, float, str)

# VCFv4.3
class varRecord:
    """This class is for storing variant information from VCF.    
    It follows VCFv4.2 specification. Data lines of VCF are 
    tab seperated, and each column represents each information
    related to variant.  For example,
    
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

    varRecord instance saves content of 'FILTER' field as list,
    content of 'INFO' field as dictionary, and content of 'FORMAT' field
    as dictionary. Other fileds are saved as string or integer or float.
    
    FILTER field : '.' -> ['.']
    INFO field   : 'AC=2;AN=2;DB;DP=182;H2;NS=65' -> {'AC': '2', 'AN': '2', ... , 'DB': '', 'NS': '65'}
    
    Each tag in INFO field has the form of 'key=value'. The explanation
    for key is at meta information lines of VCF, located at the top of VCF.
    The type of value of tag of INFO field could be 'Integer', 'Float',
    'Flag', 'Character', 'String'.

    In the case of the absence of 'FORMAT' field (and following 'SAMPLE' fields),
    'format_headers' and 'sample_info' attributes are saved as None.
    
    (optional)
    FORMAT field  : 'GT:GQ:DP:HQ' -> {1: 'GT', 2: 'GQ', 3: 'DP', 4: 'HQ'}
    
    SAMPLE fields :
    # ... FORMAT        NA00001         NA00002
    ..    GT:GQ:DP:HQ   0|0:48:1:51,51  1|0:48:8:51,51
    -> {0: 'GT', 1: 'GQ', 2: 'DP', 3: 'HQ'}
    -> {'NA00001' : {'DP': 1, 'GQ': 48, 'GT': '0|0', 'HQ': '51,51'},
        'NA00002' : {'DP': 8, 'GQ': 48, 'GT': '1|0', 'HQ': '51,51'}}
    
    Each tag in FORMAT field is concatenated by ':'. Following SAMPLE fields
    also have the form of ':', values for tags in FORMAT field are concatenated
    by ':'. The explanation for tags of FORMAT field is at meta information lines.
    The type of value of tag of FORMAT filed could be 'Integer', 'Float',
    'Character', 'String' (same as INFO field except 'Flag').
    
    Example
    -------    
    >>> sample_list = ["NA00001", "NA00002"]
    >>> test_str = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51"
    >>> test_info = test_str.split()
    >>> test_var = varRecord(*test_info[:8], format_header=test_info[8], format_list=test_info[9:], sample_name_list=sample_list)
    >>> print(test_var.chrom)
    20
    >>> print(test_var.pos)
    14370
    >>> print(test_var.ID)
    rs6054257
    >>> print(test_var.ref)
    G
    >>> print(test_var.alt)
    A
    >>> print(test_var.qual)
    29.0
    >>> print(test_var.filter)
    ['PASS']
    >>> print([ (k, test_var.info[k]) for k in sorted(list(test_var.info.keys())) ])
    [('AF', 0.5), ('DB', True), ('DP', 14), ('H2', True), ('NS', 3)]
    >>> print([ (k, test_var.format_headers[k]) for k in sorted(list(test_var.format_headers.keys())) ])
    [(0, 'GT'), (1, 'GQ'), (2, 'DP'), (3, 'HQ')]
    >>> for sample in sorted(list(test_var.sample_info.keys())):
    ...     for key in sorted(list(test_var.sample_info[sample].keys())):
    ...         print(sample, key, test_var.sample_info[sample][key])
    NA00001 DP 1
    NA00001 GQ 48
    NA00001 GT 0|0
    NA00001 HQ 51,51
    NA00002 DP 8
    NA00002 GQ 48
    NA00002 GT 1|0
    NA00002 HQ 51,51
    >>> print(test_var.get_column_num())
    10
    """
    
    __slots__ = ('_chrom', '_pos', '_ID',
                 '_ref', '_alt', '_qual',
                 '_filter', '_info',
                 '_format_headers', 'sample_info')
    
    def __init__(self,
        chrom: str,
        pos: str,
        ID: str,
        ref: str,
        alt: str,
        qual: str,
        var_filter: str,
        info: str,
        meta_info: Dict[str, dict],
        format_header: str = False,
        format_list: List[str] = False,
        sample_name_list: List[str] = False,
        ) -> None:
        
        """Init instance for each variant

        Parameters
        ----------
        chrom : str
            Chromosome
        pos : int
            Reference position (1-base coordinate)
        ID : str
            Variant identifier
        ref : str
            Reference allele
        alt : str
            Alternative allele
        qual : int
            Variant quality
        var_filter : str
            Filter status
        info : str
            Additional info combined across all samples
        meta_info : Dict[str, dict]
            Dictionary containing meta information lines
        format_header : str, optional
            Header for additional info of each sample, by default False
        format_list : list, optional
            Additional info of each sample_, by default False
        sample_name_list : list, optional
            Sample name list_, by default False
        """
        
        self.chrom = chrom
        self.pos = pos
        self.ID = ID
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = var_filter
        self._info : Dict[str, INFO] = self._parse_info(info, meta_info)
        
        ## optional attributes
        if format_header and format_list and sample_name_list:
            self.format_headers : Dict[str, str] = format_header
            self.sample_info : Dict[str, Dict[str, str]] = self._parse_sample_format(format_list, sample_name_list)
        else:
            self.format_headers = None
            self.sample_info = None

        return
    
    @property
    def chrom(self) -> str:        
        return self._chrom
    
    @chrom.setter
    def chrom(self, _chrom : str) -> None:
        self._chrom = _chrom.strip()
        return
    
    @property
    def pos(self) -> int:
        return self._pos
    
    @pos.setter
    def pos(self, _pos : str) -> None:
        try:
            self._pos = int(_pos)
        except ValueError:
            sys.exit('POS data must be int type')
        return
    
    @property
    def ID(self) -> str:
        return self._ID
    
    @ID.setter
    def ID(self, _ID : str) -> None:
        self._ID = _ID
        return
    
    @property
    def ref(self) -> str:
        return self._ref
    
    @ref.setter
    def ref(self, _ref : str) -> None:
        self._ref = _ref
        return
    
    @property
    def alt(self) -> str:
        return self._alt
    
    @alt.setter
    def alt(self, _alt : str) -> None:
        self._alt = _alt
        return

    @property
    def qual(self) -> float:
        return self._qual
    
    @qual.setter
    def qual(self, _qual : str) -> None:
        try:
            self._qual = float(_qual)
        except ValueError:
            sys.exit('QUAL data must be float type')
    
    @property
    def filter(self) -> str:
        return self._filter
    
    @filter.setter
    def filter(self, _filter : str) -> None:
        self._filter = []
        if _filter.endswith(';'):
            _filter = _filter[:-1]
        self._filter = _filter.split(';')
        
        return

    @property
    def info(self) -> Dict:
        return self._info
    
    @property
    def format_headers(self) -> Dict[str, str]:
        return self._format_headers
    
    @format_headers.setter
    def format_headers(self, _format_header : str) -> None:
        if _format_header:
            format_headers = {}
            for i, key in enumerate(_format_header.split(':')):
                format_headers[i] = key
            self._format_headers = format_headers
        else:
            self._format_headers = dict()
        return
    
    def _transform_type(self,
        value : str,
        type : str
        ) -> INFO or FORMAT:
        
        """Transform input value from string type
        to proper python data type according to the meta
        information.
        
        -----------------------------
        Types     | Python data type
        -----------------------------
        Integer   | int
        Float     | float
        Flag      | bool
        Character | str
        String    | str
        -----------------------------
        """
        
        try:
            if type == 'Integer':
                return int(value)
            elif type == 'Float':
                return float(value)
            elif type == 'Flag':
                return bool(value)
        except ValueError:
            return value        

    def _parse_info(self,
        _info : str,
        meta_info : Dict[str, dict]
        ) -> Dict[str, INFO]:
        
        """Parse 'key=value' pairs seperated by ';' of INFO field
        to dictionary. Values are stored by proper data type, according to
        the type of the key(tag).
        
        -----------------------------
        Types     | Python data type
        -----------------------------
        Integer   | int
        Float     | float
        Flag      | bool
        Character | str
        String    | str
        -----------------------------
        """
        
        self_info : Dict[str, INFO] = {}
        if _info.endswith(';'):
            _info = _info[:-1]
        for info_data in _info.split(';'):
            try:
                key, value = info_data.split('=')
            except ValueError:
                ## case that info_data is 'FLAG' type
                key = info_data
                ## if it has info_data which of type is 'FLAG',
                ## its value should be 'O' (otherwise '')
                value = 'O'
            
            self_info[key.strip()] = self._transform_type(value.strip(), meta_info['INFO'][key.strip()].type)
        
        return self_info

    def _parse_sample_format(self,
        sample_formats : List[str] = False,
        sample_names : List[str] = False
        ) -> Dict[str, Dict[str, str]]:
        
        if sample_formats:
            sample_info = { sample_name : {} for sample_name in sample_names }
            if len(sample_formats) != len(sample_names):
                sys.exit('The number of FORMAT columns and the number of samples must be the same.')
            
            for sample_i, sample_format in enumerate(sample_formats):
                sample_name = sample_names[sample_i]
                for i, value in enumerate(sample_format.split(':')):
                    sample_info[sample_name][self.format_headers[i]] = value
        return sample_info

    def __str__(self) -> str:
        res_str = f'varRecord'
        res_str += f'({self.chrom}, {self.pos}'
        res_str += f', {self.ID}, {self.ref}'
        res_str += f', {self.alt}, {self.qual}'
        res_str += f', {self.filter}, {self.info}'
        if not self.format_headers and not self.sample_info:
            res_str += ')'
            return res_str
        else:
            res_str += f', {self.format_headers}, {self.sample_info})'
            return res_str
    
    def get_column_num(self) -> int:
        """Return # of columns at VCF data line"""
        cols = [ getattr(self, c) for c in self.__slots__ if getattr(self, c) ]
        return len(cols)


def _test():
    """Test function of this code"""
    
    import doctest
    doctest.testmod()
    
    return

if __name__ == "__main__":
    _test()