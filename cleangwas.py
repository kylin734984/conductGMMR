import pandas as pd
import sys
from decimal import Decimal,getcontext
import time
import argparse
from functools import wraps
import logging
import os



'''
Summary:
    This a script that conducts basic and advanced quality control for a gwas summary data.
Usage:
    1) Configure basic variables in the head of this script.
    2) using following command: python3 cleangwas.py -h for help.
Details:
    1) SNPs meeting following conditions will be removed:beta<0
        se<0
        P<0
        info<0
        freq<0 or freq>1
        or < 0

    2) SNPs meeting will be removed if you set additional parameters:
        p>threshold
        info<threshold
        freq<maf
        N < threshold
        duplicated SNPs
        palindrome SNPs
        indels
        SNPs with missing values
    3) Alleles of SNPs will be uppercased.
    4) Uppercase RS number
    5) Eligible SNPs will be written into file you specified.    

'''

getcontext().prec = 1000
class options():
    def __init__(self):
        self.infp = 'HDL'
        self.outfp = 'HDL_test'
        self.pThresh = float('0.5')
        self.identifier = 'RS:A1:A2'
        self.maf = 0.05
        self.info = 0.8
        self.sep = '\s+'
        self.n = 0.0
        self.dropna = False
        self.rmpali = True
        self.rmindel = True
        self.truncate = 'basicgwas'
        self.unique = True
        self.include = None
        self.exclude = None
        self.rs_name = None
        self.chr_name = None
        self.pos_name = None
        self.a1_name = None
        self.a2_name = 'Aa2'
        self.freq_name = 'FREQ.A1.1000G.EUR'
        self.beta_name = None
        self.or_name = None
        self.se_name = None
        self.info_name = None
        self.P_name = None
        self.N_name = None


pd.set_option('precision', 5)
bad_snps = {
    'P0': pd.DataFrame(),  # SNPs with p<0 or p=1
    'INFO0': pd.DataFrame(),  # imputaion <0 or >1
    'FREQ0': pd.DataFrame(),  # freq <0 or >1
    'SE0': pd.DataFrame(),  # se <0
    'OR0': pd.DataFrame(),  # odds ratio <1
    'BETA0': pd.DataFrame(),  # beta<0
    'N0': pd.DataFrame(),  # n<0
    'P1': pd.DataFrame(),  # p< p threshold
    'INFO1': pd.DataFrame(),  # imputation < info threshold
    'FREQ1': pd.DataFrame(),  # frq > maf threshold
    'REPEAT0': pd.DataFrame(),  # duplicated SNPs
    'PALI0': pd.DataFrame(),  # Palindrome SNPs
    'NA': pd.DataFrame(),  # SNPs with NA
    'DUP': pd.DataFrame()  # duplicated SNPs
}


def summarise(opts, del_no, cnames, logger):
    if 'P' in cnames:
        logger.logger.info(
            '{0} SNPs with P value greater than {1} or less than 0 were deleted.'.format(del_no['P0'], opts.pThresh))
    if 'FREQ' in cnames:
        logger.logger.info('{0} SNPs with maf > {1} were deleted.'.format(del_no['FREQ0'], opts.maf))
    if 'OR' in cnames:
        logger.logger.info('{0} SNPs with OR < 0 were deleted.'.format(del_no['OR0']))
    if 'SE' in cnames:
        logger.logger.info('{0} SNPs with SE < 0 were deleted.'.format(del_no['SE0']))
    if 'N' in cnames:
        logger.logger.info('{0} SNPs with sample size < {1} were deleted.'.format(del_no['N0'], opts.n))
    if 'INFO' in cnames:
        logger.logger.info('{0} SNPs with imputation score < {1} or with missing value were deleted.'.format(del_no['INFO0'], str(opts.info)))
    if opts.dropna:
        logger.logger.info('{0} SNP with missing value were deleted.'.format(del_no['NA']))
    if opts.rmindel:
        logger.logger.info('{0} indel SNP were deleted.'.format(del_no['INDEL']))
    if opts.rmpali:
        if 'A1' in cnames and 'A2' in cnames:
            logger.logger.info('{0} palindromic SNPs were deleted.'.format(del_no['PALI0']))
    if opts.unique:
        if 'A1:A2' in opts.identifier and not opts.rmpali:
            logger.logger.info('{0} palindromic SNPs were deleted.'.format(del_no['PALI0']))
        logger.logger.info('{0} duplicated SNPs were deleted.'.format(del_no['DUP']))


del_no = {
    'P0': 0,  # SNPs with p<0 or p=1
    'INFO0': 0,  # imputaion <0 or < threshold
    'FREQ0': 0,  # freq <0 or >1
    'SE0': 0,  # se <0
    'OR0': 0,  # odds ratio <1
    'BETA0': 0,  # beta<0
    'N0': 0,  # n<0
    'P1': 0,  # p< p threshold
    'INFO1': 0,  # imputation < info threshold
    'FREQ1': 0,  # frq > maf threshold
    'PALI0': 0,  # Palindrome SNPs
    'NA': 0,  # SNPs with NA
    'DUP': 0,  # duplicated SNPs
    'INDEL': 0
}
default_cnames = {

    # RS NUMBER
    'SNP': 'RS',
    'MARKERNAME': 'RS',
    'SNPID': 'RS',
    'RS': 'RS',
    'RSID': 'RS',
    'RS_NUMBER': 'RS',
    'RS_NUMBERS': 'RS',
    #CHR
    'CHR':'CHR',
    'CHROMOSOME':'CHR',
    #POS
    'POS':'POS',
    'BP':'POS',
    'POSITION':'POS',
    # NUMBER OF STUDIES
    'NSTUDY': 'NSTUDY',
    'N_STUDY': 'NSTUDY',
    'NSTUDIES': 'NSTUDY',
    'N_STUDIES': 'NSTUDY',
    # P-VALUE
    'P': 'P',
    'PVALUE': 'P',
    'P_VALUE': 'P',
    'PVAL': 'P',
    'P_VAL': 'P',
    'GC_PVALUE': 'P',
    # ALLELE 1
    'A1': 'A1',
    'ALLELE1': 'A1',
    'ALLELE_1': 'A1',
    'EFFECT_ALLELE': 'A1',
    'REFERENCE_ALLELE': 'A1',
    'INC_ALLELE': 'A1',
    'EA': 'A1',
    # ALLELE 2
    'A2': 'A2',
    'ALLELE2': 'A2',
    'ALLELE_2': 'A2',
    'OTHER_ALLELE': 'A2',
    'NON_EFFECT_ALLELE': 'A2',
    'DEC_ALLELE': 'A2',
    'NEA': 'A2',
    # N
    'N': 'N',
    'NCASE': 'N_CAS',
    'CASES_N': 'N_CAS',
    'N_CASE': 'N_CAS',
    'N_CASES': 'N_CAS',
    'N_CONTROLS': 'N_CON',
    'N_CAS': 'N_CAS',
    'N_CON': 'N_CON',
    'N_CASE': 'N_CAS',
    'NCONTROL': 'N_CON',
    'CONTROLS_N': 'N_CON',
    'N_CONTROL': 'N_CON',
    'WEIGHT': 'N',  # metal does this. possibly risky.
    # SIGNED STATISTICS
    'ZSCORE': 'Z',
    'Z-SCORE': 'Z',
    'GC_ZSCORE': 'Z',
    'Z': 'Z',
    'Z_VALUE':'Z',
    'ZVALUE':'Z',
    'ZVAL':'Z',
    'OR': 'OR',
    'B': 'BETA',
    'BETA': 'BETA',
    'LOG_ODDS': 'LOG_ODDS',
    'EFFECTS': 'BETA',
    'EFFECT': 'BETA',
    'SIGNED_SUMSTAT': 'SIGNED_SUMSTAT',
    # INFO
    'INFO': 'INFO',
    # MAF
    'EAF': 'FREQ',
    'FRQ': 'FREQ',
    'MAF': 'FREQ',
    'FRQ_U': 'FREQ',
    'FREQ': 'FREQ',
    'FREQ1': 'FREQ',
    # se
    'SE': 'SE',
    'STDERR': 'SE'
}


class Logger():
    def __init__(self, name=None):
        self.logger = logging.getLogger(name) if name else logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.formatter = logging.Formatter('%(message)s')
        self.f_handler = logging.FileHandler('{0}.log'.format(__name__))
        self.f_handler.setFormatter(self.formatter)
        self.f_handler.setLevel(logging.INFO)
        self.logger.addHandler(self.f_handler)
        self.console_handler = logging.StreamHandler()
        self.console_handler.setFormatter(self.formatter)
        self.console_handler.setLevel(logging.INFO)
        self.logger.addHandler(self.console_handler)
        self.logger.setLevel(logging.INFO)

    def rename(self, outprefix):
        self.f_handler.flush()
        self.f_handler.close()
        self.logger.removeHandler(self.f_handler)
        if isinstance(outprefix, str):
            try:
                os.remove('{0}.log'.format(outprefix))
            except:
                pass
            try:
                os.rename('{0}.log'.format(__name__), '{0}.log'.format(outprefix))
            except:
                print('Can\'t rename log file: {0}'.format('{0}.log'.format(__name__)))


logger = Logger()

def timeit(logger):
    def decorator(fun):
        @wraps(fun)  # 为了保留被装饰函数的函数名和帮助文档信息
        def wrapper(*args, **kwargs):
            """这是一个wrapper函数"""
            start_time = time.time()
            logger.logger.info('Start time: {0}'.format(time.strftime('%H:%M:%S', time.localtime(start_time))))
            res = fun(*args, **kwargs)
            end_time = time.time()
            use = time.gmtime(end_time - start_time)
            logger.logger.info(
                'Complete time: {0}\nEclipsed time: {1}'.format(time.strftime('%H:%M:%S', time.localtime(end_time)),
                                                                time.strftime('%H:%M:%S', use)))
            return res
        return wrapper
    return decorator

def parse_arg():
    parser = argparse.ArgumentParser(description='Clean GWAS data.')
    parser.add_argument('-i', '--input', dest='infp', required=True, default=None, metavar='FILE', help='Input file')
    parser.add_argument('-o', '--output', dest='outfp', required=True, default=None, metavar='FILE', help='Output file')
    parser.add_argument('--pThresh', default=1, type=float, metavar='0-1',
                        help='Maximum P value threshold for containing')
    parser.add_argument('--maf', default=0, type=float, metavar='0-0.5',
                        help='Minor allele frequency threshold for retain')
    parser.add_argument('--info', default=0, type=float, metavar='0-1', help='imputation score threshold')
    parser.add_argument('--n', default=0, type=float, metavar='>0', help='Specify a minimum sample size to restore')
    parser.add_argument('--sep', default='\s+', metavar='separator', help='delimiter of input file')
    parser.add_argument('--dropna', action='store_true', help='Drop records with NA values')
    parser.add_argument('--rmpalindromic', dest= 'rmpali', action='store_true', help='Drop palindromic SNPs')
    parser.add_argument('--rmindel', action='store_true', help='Drop indel SNPs')
    parser.add_argument('--truncate', choices=['basicgwas'], help='Create a subset.'
                                                                  '\n[basicgwas] represents create a gwas inculding RS A1 A2 FREQ BETA(OR) SE P N')
    group = parser.add_argument_group()
    group.add_argument('--unique', action='store_true', help='Drop duplicated SNPs depending on --identifier parameter')
    group.add_argument('--identifier', default='RS',
                       choices=['RS:A1:A2', 'RS:CHR:POS', 'RS', 'RS:CHR:POS:A1:A2', 'CHR:POS:A1:A2'],
                       help='SNP identity used for denote a unique marker')
    parser.add_argument('--include', default=None, metavar='FILE', help='SNP list file without header includes RS# you need. ')
    parser.add_argument('--exclude', default=None, metavar='FILE', help='SNP list file without header includes RS# you did not need')
    parser.add_argument('--rs-name', help='Specify rs column name')
    parser.add_argument('--chr-name', help='Specify chromosome column name')
    parser.add_argument('--pos-name', help='Specify position column name')
    parser.add_argument('--a1-name', help='Specify effect allele column name')
    parser.add_argument('--a2-name', help='Specify non effect allele column name')
    parser.add_argument('--freq-name', help='Specify frequency column name')
    parser.add_argument('--beta-name', help='Specify beta coefficient column name')
    parser.add_argument('--or-name', help='Specify OR coefficient column name')
    parser.add_argument('--se-name', help='Specify SD column name')
    parser.add_argument('--info-name', help='Specify imputation score column name')
    parser.add_argument('--N-name', help='Specify sample size column name')
    parser.add_argument('--P-name', help='Specify P value column name')

    return parser.parse_args()


def read_header(f_p: str, sep='\s+') -> str:
    df = pd.read_csv(f_p, sep=sep, header=0, nrows=1)
    return '\t'.join(list(df.columns))


def parse_header(header: str, opts, logger, default_cnames=default_cnames) -> list:
    namedict = {
        opts.rs_name: 'RS',
        opts.chr_name: 'CHR',
        opts.pos_name: 'POS',
        opts.a1_name: 'A1',
        opts.a2_name: 'A2',
        opts.freq_name: 'FREQ',
        opts.beta_name: 'BETA',
        opts.or_name: 'OR',
        opts.se_name: 'SE',
        opts.info_name: 'INFO',
        opts.N_name: 'N',
        opts.P_name: 'P'
    }
    namedict.pop(None, '1')
    namedict = {k.replace('-', '_').replace('.', '_').upper(): v for k, v in namedict.items()}
    default_cnames.update(namedict)
    he_sep = header.split()
    new_he_sep = []
    notes = 'Parsing header for {0}:\n'.format(opts.infp)
    for e in he_sep:
        ewrapper = e.upper().replace('-', '_').replace('.', '_')
        if ewrapper in default_cnames.keys():
            new_he_sep.append(default_cnames.get(ewrapper))
            notes += 'Identified:   {0}{1}{2}\n'.format(e, ' ' * (20 - len(e)), default_cnames.get(ewrapper))
        else:
            new_he_sep.append(e)
            notes += 'unidentified: {0}{1}Unknown\n'.format(e, ' ' * (20 - len(e)))
    logger.logger.info(notes)
    logger.logger.info('Tips: If you use some unidentified columns, you will get an error!')
    return new_he_sep


def split_df(df: pd.DataFrame, flag: pd.Series) -> tuple:
    return df.loc[flag].copy(), df.loc[flag == False].copy()


def fliter_p(a_row: pd.Series, p_thresh=1) -> bool:
    return False if Decimal.is_nan(a_row['P']) else 0 <= Decimal(a_row['P']) <= p_thresh


def fliter_info(a_row: pd.Series, info_thresh=0) -> bool:
    return False if Decimal.is_nan(Decimal(a_row['INFO'])) else info_thresh <= Decimal(a_row['INFO']) <= 1


def fliter_freq(a_row: pd.Series, maf_thresh=0) -> bool:
    return False if Decimal.is_nan(a_row['FREQ']) else maf_thresh <= Decimal(a_row['FREQ']) <= 1 - maf_thresh


def fliter_se(a_row: pd.Series) -> bool:
    return False if Decimal.is_nan(a_row['SE']) else Decimal(a_row['SE']) > 0

def fliter_or(a_row: pd.Series) -> bool:
    return False if Decimal.is_nan(a_row['OR']) else Decimal(a_row['OR']) > 0

def fliter_n(a_row: pd.Series, N=0) -> bool:
    return False if Decimal.is_nan(a_row['N']) else Decimal(a_row['N']) > N


def rmindel(a_row: pd.Series) -> bool:
    return len(a_row['A1'] + a_row['A2']) > 2


def is_palindromic(arow):
    return (arow['A1'] + ':' + arow['A2']).upper() in ['A:T', 'T:A', 'C:G', 'G:C']

def is_dup(a_row, idpool, identifier):
    id_ = ':'.join([a_row[x] for x in identifier.strip().split(":")])
    return True if idpool.isin([id_]).tolist().count(True) > 1 else False

def get_converter(cnames):
    '''
    CHR 因为可能有XY染色体，所以当做str类型
    '''
    converter = {}
    dec_cols = ['FREQ','BETA', 'SE','P', 'N', 'OR', 'Z', 'POS']
    for x in cnames:
        converter[x] = Decimal if x in dec_cols else str
    return converter


def split_dup(df, opts, cnames):
    usednames = opts.identifier.strip().split(':')
    if set(usednames) - set(cnames):
        logger.logger.warning(
            '\033[5;31mWarning: --indentifier use columns unidentified in your file:{0}, Please check it.\nExit...  \033[0m'.format(
                set(usednames) - set(cnames)))
        sys.exit(1)
    id_prefix = opts.identifier.strip().replace('A1:A2', '').strip(':').split(":")
    dup_df = df.loc[df.duplicated(subset=id_prefix, keep=False)].copy()
    if dup_df.empty:
        return df, pd.DataFrame()
    idpool1 = dup_df.apply(lambda x: ':'.join([x[e] for e in usednames]).upper(), axis=1)
    if 'A1:A2' in opts.identifier:
        usednames = opts.identifier.replace('1', '3').replace('2', '1').replace('3', '2').split(':')
        idpool2 = dup_df.apply(lambda x: ':'.join([x[e] for e in usednames]).upper(), axis=1)
        cnames2 = opts.identifier.replace('1', '3').replace('2', '1').replace('3', '2').split(':')
        # dup_df = dup_df.loc[:, cnames2].copy()
        # idpool2 = dup_df.apply(lambda x: ':'.join([e for e in x]).upper(), axis=1)
        idpool3 = idpool1.apply(
            lambda x: x.replace('A', 'Z').replace('T', 'A').replace('Z', 'T').replace('C', 'Z').replace('G',
                                                                                                        'C').replace(
                'Z', 'G'))
        idpool4 = idpool2.apply(
            lambda x: x.replace('A', 'Z').replace('T', 'A').replace('Z', 'T').replace('C', 'Z').replace('G',
                                                                                                        'C').replace(
                'Z', 'G'))
        idpool1 = pd.concat([idpool1, idpool2, idpool3, idpool4], ignore_index=True, sort=False)
    drop_dup_df, remain_dup_df = split_df(dup_df,
                                          dup_df.apply(is_dup, axis=1, idpool=idpool1, identifier=opts.identifier))
    df = df.append(drop_dup_df, ignore_index=True, sort=False).drop_duplicates(keep=False)
    return df, drop_dup_df


def qc(opts, cnames: list, logger=logger) -> pd.DataFrame:
    total_df = pd.DataFrame()
    converter = get_converter(cnames)
    for chunk in pd.read_csv(opts.infp, sep=opts.sep, header=0, na_values='NA', names=cnames, dtype=str,  iterator=True, chunksize=1000000):
        '''不在read_csv()使用converters参数是因为当文件最后一列是空值时，会报decimal.InvalidOperation: [<class 'decimal.ConversionSyntax'>]。
        另外同时用names和converters选项会报警告。'''
        for k,v in converter.items(): chunk[k] = chunk[k].apply(v)
        try:
            if 'P' in chunk.columns:
                chunk, dropdf = split_df(chunk, chunk.apply(fliter_p, axis=1, p_thresh=opts.pThresh))
                del_no['P0'] += dropdf.shape[0]
            elif not opts.pThresh == 1:
                logger.logger.warning(
                    '\033[5;31mWarning: --pThresh is not valid, as P column is not detected in your file \033[0m ')
            if 'RS' in chunk.columns: chunk['RS'] = chunk['RS'].str.upper()
            if 'A1' in chunk.columns: chunk['A1'] = chunk['A1'].str.upper()
            if 'A2' in chunk.columns: chunk['A2'] = chunk['A2'].str.upper()
            if opts.rmindel:
                if 'A1' in chunk.columns and 'A2' in chunk.columns:
                    dropdf, chunk = split_df(chunk, chunk.apply(rmindel, axis=1))
                    del_no['INDEL'] += dropdf.shape[0]
                else:
                    logger.logger.warning(
                        '\033[5;31mWarning: Allele columns are not found in your file, can\'t remove indel SNP  \033[0m')
            if 'FREQ' in chunk.columns:
                chunk, dropdf = split_df(chunk, chunk.apply(fliter_freq, axis=1, maf_thresh=Decimal.from_float(opts.maf)))
                # drop_snps['FREQ0'] = drop_snps['FREQ0'].append(dropdf)
                del_no['FREQ0'] += dropdf.shape[0]
            if 'OR' in chunk.columns:
                chunk, dropdf = split_df(chunk, chunk.apply(fliter_or, axis=1))
                # drop_snps['OR0'] = drop_snps['OR0'].append(dropdf)
                del_no['OR0'] += dropdf.shape[0]
            if 'SE' in chunk.columns:
                chunk, dropdf = split_df(chunk, chunk.apply(fliter_se, axis=1))
                # drop_snps['SE0'] = drop_snps['SE0'].append(dropdf)
                del_no['SE0'] += dropdf.shape[0]
            if 'N' in chunk.columns:
                chunk, dropdf = split_df(chunk, chunk.apply(fliter_n, axis=1, N=opts.n))
                # drop_snps['N0'] = drop_snps['N0'].append(dropdf)
                del_no['N0'] += dropdf.shape[0]
            if 'INFO' in chunk.columns:
                chunk, dropdf = split_df(chunk, chunk.apply(fliter_info, axis=1, info_thresh=opts.info))
                # drop_snps['INFO0'] = drop_snps['INFO0'].append(dropdf)
                del_no['INFO0'] += dropdf.shape[0]
            total_df = total_df.append(chunk, sort=False)
        except ValueError:
            if not chunk.empty:
                logger.logger.info('Check qc code!!, something may be wrong')
    if opts.rmpali and not total_df.empty:
        if 'A1' in total_df.columns and 'A2' in total_df.columns:
            dropdf, total_df = split_df(total_df, total_df.apply(is_palindromic, axis=1))
            del_no['PALI0'] += dropdf.shape[0]
        else:
            logger.logger.warning(
                '\033[5;31mWarning: Can not delete palindromic SNPs because of missing alleles columns\nExit...\033[0m')
            sys.exit(1)
    elif total_df.empty:
        logger.logger.warning('Warning: no eligible SNPs left.')
    if opts.unique and not total_df.empty:
        if 'A1:A2' in opts.identifier and not opts.rmpali:
            dropdf, total_df = split_df(total_df, total_df.apply(is_palindromic, axis=1))
            del_no['PALI0'] += dropdf.shape[0]
            logger.logger.info(
                'As you specified --unique option and use Allele columns, palindromic SNPs will be deleted regardless of whether specifing --rmpali .')
        total_df, dropdf = split_dup(total_df, opts, cnames)
        del_no['DUP'] += dropdf.shape[0]
    elif total_df.empty:
        logger.logger.info('Warning: no eligible SNPs left.')

    summarise(opts, del_no, cnames, logger)
    return total_df


def resort_col(df: pd.DataFrame) -> pd.DataFrame:
    '''
    将一个dataframe重新排序。'RS', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N'在前，其他在后。
    '''
    std_cols = ['RS', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N', 'OR']
    act_cols = [e for e in std_cols if e in df.columns]
    add_cols = [e for e in df.columns if not e in std_cols]
    act_cols.extend(add_cols)
    df = df[act_cols].copy()
    return df


def get_identifier_pool(opts, headercols, logger):
    usednames = opts.identifier.split(':')
    if set(usednames) - set(headercols):
        logger.logger.warning(
            '\033[5;31mWarning: --indentifier use columns unidentified in your file:{0}, Please check it.\nExit...  \033[0m'.format(
                set(usednames) - set(headercols)))
        sys.exit(1)
    df = pd.read_csv(opts.infp, names=headercols, usecols=usednames, sep=opts.sep, header=0)
    bad_snps['NA'] = bad_snps['NA'].append(df[df.isna().apply(lambda x: any(x), axis=1)], sort=False)
    df.dropna(inplace=True)
    idpool1 = df.apply(lambda x: ':'.join([e for e in x]).upper(), axis=1)
    if 'A1:A2' in opts.identifier:
        cnames2 = opts.identifier.replace('1', '3').replace('2', '1').replace('3', '2').split(':')
        df = df.loc[:, cnames2].copy()
        idpool2 = df.apply(lambda x: ':'.join([e for e in x]).upper(), axis=1)
        idpool3 = idpool1.apply(
            lambda x: x.replace('A', 'Z').replace('T', 'A').replace('Z', 'T').replace('C', 'Z').replace('G',
                                                                                                        'C').replace(
                'Z', 'G'))
        idpool4 = idpool2.apply(
            lambda x: x.replace('A', 'Z').replace('T', 'A').replace('Z', 'T').replace('C', 'Z').replace('G',
                                                                                                        'C').replace(
                'Z', 'G'))
        idpool1 = pd.concat([idpool1, idpool2, idpool3, idpool4], ignore_index=True)
    return idpool1


def dropna(opts, df: pd.DataFrame) -> pd.DataFrame:
    if opts.dropna:
        dropdf, df = split_df(df, df.isna().apply(lambda x: any(x), axis=1))
        del_no['NA'] += dropdf.shape[0]
    return df


def truncate(opts, df: pd.DataFrame) -> pd.DataFrame:
    if opts.truncate == 'basicgwas':
        std_cols = ['RS', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N']
        act_cols = [e for e in std_cols if e in df.columns]
        df = df[act_cols].copy()
    return df


def selectSNP(opts, df:pd.DataFrame, logger) -> pd.DataFrame:
    if opts.include:
        if os.path.exists(opts.include):
            insnps = pd.read_csv(opts.include, header=None, names=["RS"], sep='\s+')
            insnps['RS'] = insnps['RS'].str.upper()
            df = df[df['RS'].isin(insnps['RS'])]
            if df.empty:
                logger.logger.warning('Warning: No common SNPs\nExiting...')
                sys.exit(1)
        else:
            logger.logger.warning("Warning: File: {0} does not exist. Can not excute --include parameter!".format(opts.include))
    if opts.exclude:
        if os.path.exists(opts.exclude):
            exsnps = pd.read_csv(opts.exclude, header=None, names=['RS'], sep='\s+')
            exsnps['RS'] = exsnps['RS'].str.upper()
            df = df[df['RS'].isin(exsnps['RS'])==False]
        else:
            logger.logger.warning("Warning: File: {0} does not exist. Can not excute --exclude parameter!".format(opts.exclude))

    return df


def write_file(df, opts):
    df.to_csv(opts.outfp, sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)


@timeit(logger)
def init(opts):
    header = read_header(opts.infp, opts.sep)
    cnames = parse_header(header, opts, logger, default_cnames)
    df = qc(opts, cnames, logger)
    df = selectSNP(opts, df, logger)
    df = resort_col(df)
    df = truncate(opts, df)
    write_file(df, opts)



if __name__ == '__main__':
    if len(sys.argv) == 1:
        opts = options()
    else:
        opts = parse_arg()
    init(opts)
    logger.logger.info('End!')
    logger.rename(opts.outfp)
    sys.exit(0)
