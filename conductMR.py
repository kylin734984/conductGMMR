import argparse
import logging
import sys
import os
from functools import wraps,reduce
import time
import cleangwas
import copy
import pandas as pd
import tempfile
import subprocess
from io import StringIO
from platform import system as system_platform
#import rpy2.robjects as robjects

class options():
    def __int__(self):
        self.infp = None
        self.outfp = None

class Logger():
    def __init__(self, name = None):
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
        if isinstance(outprefix,str):
            try:
                os.remove('{0}.log'.format(outprefix))
            except:
                pass
            try:
                os.rename('{0}.log'.format(__name__), '{0}.log'.format(outprefix))
            except:
                print('Can\'t rename log file: {0}'.format('{0}.log'.format(__name__)))



def timeit(logger):
    def decorator(fun):
        @wraps(fun)  #为了保留被装饰函数的函数名和帮助文档信息
        def wrapper(*args,**kwargs):
            """这是一个wrapper函数"""
            start_time = time.time()
            logger.logger.info('Start time: {0}'.format(time.strftime('%H:%M:%S', time.localtime(start_time))))
            res = fun(*args,**kwargs)
            end_time = time.time()
            use = time.gmtime(end_time - start_time)
            logger.logger.info('Complete time: {0}\nEclipsed time: {1}'.format(time.strftime('%H:%M:%S', time.localtime(end_time)),
                                                                  time.strftime('%H:%M:%S', use)))
            return res
        return wrapper
    return decorator

def parse_arg():
    parser = argparse.ArgumentParser(description='A module API demo')
    exp_group = parser.add_mutually_exclusive_group(required=True)
    exp_group.add_argument('--exposure', default=None, metavar='FILE', help='exposure file')
    exp_group.add_argument('--exposure-gm', default=None, metavar='FILE', help='exposure file for gut microbia')
    parser.add_argument('--outcome', required=True, default=None, metavar='FILE', help='outcome file')
    parser.add_argument('-o', '--out', required=True, dest='outfp', metavar="file prefix")
    parser.add_argument('--pThresh', default=1, type=float, metavar='0-1',
                        help='Maximum P value threshold of SNPs for exposure')
    parser.add_argument('--maf', default=0, type=float, metavar='0-0.5',
                        help='Minor allele frequency threshold for retain')
    parser.add_argument('--info', default=0, type=float, metavar='0-1', help='imputation score threshold')
    parser.add_argument('--n', default=0, type=float, metavar='>0', help='Specify a minimum sample size to restore')
    parser.add_argument('--sep', default='\s+', metavar='separator', help='delimiter of input file')
    parser.add_argument('--dropna', action='store_true', help='Drop records with NA values')
    parser.add_argument('--rmpalindromic', dest= 'rmpali', action='store_true', help='Drop palindromic SNPs')
    parser.add_argument('--rmindel', action='store_true', help='Drop indel SNPs')
    parser.add_argument('--truncate', default='basicgwas', choices=['basicgwas'], help='Create a subset.'
                                                                  '\n[basicgwas] represents create a gwas inculding RS A1 A2 FREQ BETA(OR) SE P N')
    group = parser.add_argument_group()
    group.add_argument('--unique', action='store_true', help='Drop duplicated SNPs depending on --identifier parameter')
    group.add_argument('--identifier', default='RS:A1:A2',
                       choices=['RS:A1:A2', 'RS:CHR:POS', 'RS', 'RS:CHR:POS:A1:A2', 'CHR:POS:A1:A2'],
                       help='SNP identity used for denote a unique marker')
    parser.add_argument('--include', default=None, metavar='FILE',
                        help='SNP list file without header includes RS# you need. ')
    parser.add_argument('--exclude', default=None, metavar='FILE',
                        help='SNP list file without header includes RS# you did not need')
    parser.add_argument('--mr-method', action='append', default=['mr_wald_ratio','mr_ivw', 'mr_egger_regression'],
                        choices=['mr_wald_ratio', 'mr_two_sample_ml','mr_egger_regression', 'mr_weighted_median', 'mr_ivw', 'mr_presso', 'mr_gsmr'], help='MR method you need.')
    parser.add_argument('--heterogeneity', default='mr_presso', choices=['mr_presso', 'mr_gsmr', 'none'], help='Pleiotropy test method ')
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
    parser.add_argument('--exp-name', default='exposure', help="Specify a phenotype name for exposure. This option is not in effect if you use --exposure-gm option at the same time.")
    parser.add_argument('--out-name', default='outcome', help='Specify a phenotype name for outcome')
    group1 = parser.add_argument_group()
    group1.add_argument('--clump', action='store_true', help='PLNK clumping algorithm')
    group1.add_argument('--windowsize', default=250, help='windowsize parameter for PLINK clumping')
    group1.add_argument('--r2', default=0.1, help='LD threshold for PLINK clumping')
    group1.add_argument('--clump-p1', default=1, help='p1 threshold for plink clumping')
    group1.add_argument('--origin', default='offline', help='genotype file directory')

    return parser.parse_args()

def read_exp(opts, logger):
    opts = copy.deepcopy(opts)
    opts.infp = opts.exposure
    header = cleangwas.read_header(opts.infp, opts.sep)
    cnames = cleangwas.parse_header(header, opts, logger, cleangwas.default_cnames)
    df = cleangwas.qc(opts, cnames, logger)
    df = cleangwas.selectSNP(opts, df, logger)
    df = cleangwas.resort_col(df)
    df = cleangwas.truncate(opts, df)
    return df

def read_outc(opts, logger):
    for k, v in cleangwas.del_no.items():cleangwas.del_no[k]=0
    opts = copy.deepcopy(opts)
    opts.infp = opts.outcome
    opts.pThresh = 1
    opts.include = 'include'
    header = cleangwas.read_header(opts.infp, opts.sep)
    cnames = cleangwas.parse_header(header, opts, logger, cleangwas.default_cnames)
    df = cleangwas.qc(opts, cnames, logger)
    df = cleangwas.selectSNP(opts, df, logger)
    df = cleangwas.resort_col(df)
    df = cleangwas.truncate(opts, df)
    try:
        os.remove('include')
    except:
        pass
    return df

def clumping(opts, expdf, refdir, logger):
    plinkreport = expdf.loc[:,['RS', 'P']]
    plinkreport['RS'] = plinkreport['RS'].str.lower()
    plinkreport.to_csv('plinkreport', sep='\t', header=['SNP', 'P'], index=False)
    independent = []
    if opts.clump and opts.origin == 'offline':
        for i in range(1, 23):
            try:
                errf = tempfile.TemporaryFile()
                sysstr = system_platform()
                plinkscript = os.path.join(refdir, 'refsource', 'plink_{0}'.format(sysstr), 'plink')
                geno = os.path.join(refdir, 'refsource', '1000G_EUR_Phase3_plink', '1000G.EUR.QC.{0}'.format(str(i)))
                content = subprocess.check_output('{0} --bfile {1} --clump {2} --clump-p1 {3} --clump-r2 {4} --clump-kb {5} '.format(plinkscript, geno, 'plinkreport', opts.clump_p1, opts.r2, opts.windowsize),
                                                  shell=True, stderr=errf, bufsize=0)
                errf.flush()
                errf.seek(0)
                if 'Warning: No significant --clump results' in errf.read().decode('utf-8'):
                    logger.logger.info('No significant SNPs in chromosome {0}'.format(str(i)))
                    continue
                else:
                    pass
                    #logger.logger.info(content.decode('utf-8')[-82:-35])
                errf.close()
            except subprocess.CalledProcessError as err1:

                logger.logger.warning('Warning: plink clumping algorithm for chromosome {0} failed. Please check .log file'.format(str(i)))
                errf.flush()
                errf.seek(0)
                logger.logger.warning(errf.read().decode('utf-8'))
                continue
            with open('plink.clumped', 'r') as reader:
                snps = [x.strip().split()[2] for x in reader.readlines()[1:] if not x.isspace()]
            independent.extend(snps)
        logger.logger.info("{0} independent SNPs after clumping".format(len(independent)))
        try:
            os.remove('plink.clumped')
            os.remove('plinkreport')
            os.remove('plink.log')
            os.remove('plink.nosex')
        except:
            pass
    elif opts.clump and opts.origin == 'online':
        logger.logger.info('Not developed now.\nExiting...')
        exit(1)
        pass
        # TODO 在线clumping
    else:
        independent = list(expdf['RS'])
    pd.Series(independent).to_csv('include', header=False,sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    return independent

def pleiotropy(opts, logger):
    r_errf = tempfile.TemporaryFile(mode='w+', encoding='utf-8')
    try:
        os.remove('mr.pleiosnp')
    except:
        pass
    if opts.heterogeneity == 'none':
        logger.logger.info("Notes: no pleiotropy test is adopted!")
    elif opts.heterogeneity == 'mr_presso':
        logger.logger.info('Pleiotropy test for {0} using {1}'.format(opts.exp_name, opts.heterogeneity))
        try:
            output = subprocess.check_output(
                'Rscript pleiotropy_test.R {0} {1} {2} {3} {4}'.format('exposure.txt', 'outcome.txt', opts.exp_name,
                                                               opts.out_name, opts.heterogeneity),
                stderr=r_errf, shell=True, bufsize=0)
            #logger.logger.info(output.decode('utf-8'))
            presso_pleiosnp = pd.read_csv('mr.pleiosnp', sep='\s+', header=0, dtype=str)
            if not presso_pleiosnp.empty:
                logger.logger.info(
                    'Pleiotropic SNPs for {0} and {1} using MR-PRESSO global test!'.format(opts.exp_name, opts.out_name))
                logger.logger.info(str(presso_pleiosnp))
            else:
                logger.logger.info(
                    'No pleiotropic SNPs for {0} and {1} identified by MR-PRESSO global test.'.format(opts.exp_name,  opts.out_name))
        except subprocess.CalledProcessError as e:
            r_errf.flush()
            r_errf.seek(0)
            logger.logger.info('Error info:\n{0}'.format(r_errf.read()))
            logger.logger.info("Pleiotropy test Failed.\nSkiping...")
            r_errf.seek(0)
            r_errf.truncate()
        except FileNotFoundError as e1:
            pass
    elif opts.heterogeneity == 'mr_gsmr':
        logger.logger.info('Pleiotropy test for {0} using {1}'.format(opts.exp_name, opts.heterogeneity))
        if pd.read_csv('exposure.txt', sep='\s+',header=0, dtype=str).shape[0] < 5:
            logger.logger.info('No enough SNPs to perform GSMR HEIDI outlier test')
            return
        refdir = os.path.split(os.path.realpath(__file__))[0]
        sysstr = system_platform()
        gctascript = os.path.join(refdir, 'refsource', 'gcta_{0}'.format(sysstr), 'gcta64')
        try:
            with open('gsmr_ref_data.txt', 'w') as writer:
                writer.write('\n'.join(
                    [os.path.join(refdir, 'refsource', '1000G_EUR_Phase3_plink', '1000G.EUR.QC.{0}'.format(str(i)))
                     for i in range(1, 23)]))
            with open('gsmr_exposure.txt', 'w') as writer:
                writer.write('\t'.join([opts.exp_name, 'exposure.txt']) + '\n')
            with open('gsmr_outcome.txt', 'w') as writer:
                writer.write('\t'.join([opts.out_name, 'outcome.txt']) + '\n')
            command = '{0} --mbfile gsmr_ref_data.txt --gsmr-file gsmr_exposure.txt gsmr_outcome.txt --gsmr-direction 0 --out gsmr_result_nonpleio --heidi-thresh 0.01  --diff-freq 1 --clump-r2 1 --gwas-thresh 1 --gsmr-snp-min 5'.format(gctascript)
            content2 = subprocess.check_output(command, stderr=r_errf, shell=True, bufsize=0)
            #logger.logger.info(content2.decode('utf-8'))
            if os.path.exists('gsmr_result_nonpleio.pleio_snps'):
                gsmr_pleiosnp = pd.read_csv('gsmr_result_nonpleio.pleio_snps', sep='\s+', header=None,
                                      names=['exposure', 'outcome', 'RS'], dtype=str)
                gsmr_pleiosnp['method'] = 'gsmr-outlier-correction'*gsmr_pleiosnp.shape[0]
                logger.logger.info('Pleiotropic SNPs for {0} and {1} using GSMR HEIDI outlier test!'.format(opts.exp_name, opts.out_name))
                gsmr_pleiosnp = gsmr_pleiosnp[['exposure', 'outcome','method', 'RS']]
                logger.logger.info(str(gsmr_pleiosnp))
                gsmr_pleiosnp.to_csv('mr.pleiosnp', header=True, index=False, sep='\t',na_rep='NA', float_format='%g', encoding='utf-8')
                for x in ['exposure.txt', 'outcome.txt']:
                    df = pd.read_csv(x, sep='\s+',header=0, dtype=str)
                    df = df[df['RS'].isin(gsmr_pleiosnp['RS'])==False]
                    df.to_csv(x, header=True, sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
            else:
                logger.logger.info('No pleiotropic SNPs for {0} and {1} identified by GSMR HEIDI outlier test'.format(opts.exp_name, opts.out_name))
        except subprocess.CalledProcessError as e:
            r_errf.flush()
            r_errf.seek(0)
            logger.logger.info('Error info:\n{0}'.format(r_errf.read()))
            logger.logger.info("Pleiotropy test Failed.\nSkiping...")
            r_errf.seek(0)
            r_errf.truncate()
        try:
            os.remove('gsmr_ref_data.txt')
            os.remove('gsmr_exposure.txt')
            os.remove('gsmr_outcome.txt')
            os.remove('gsmr_result_nonpleio.pleio_snps')
            os.remove('gsmr_result_nonpleio.gsmr')
            os.remove('gsmr_result_nonpleio.gsmr')
            os.remove('gsmr_result_nonpleio.badsnps')
        except:
            pass
        r_errf.close()

    else:
        logger.logger.info('No pleiotropy test for {0} using {1}'.format(opts.exp_name, opts.heterogeneity))
def read_gm(opts, logger):
    opts = copy.deepcopy(opts)
    opts.infp = opts.exposure_gm
    header = cleangwas.read_header(opts.infp, opts.sep)
    cnames = cleangwas.parse_header(header, opts, logger, cleangwas.default_cnames)
    opts.rmpali = False
    opts.rmindel = False
    opts.unique = False
    df = cleangwas.qc(opts, cnames, logger)
    return df

def read_outc(opts, logger):
    for k, v in cleangwas.del_no.items():cleangwas.del_no[k]=0
    opts = copy.deepcopy(opts)
    opts.infp = opts.outcome
    opts.pThresh = 1
    opts.include = 'include'
    header = cleangwas.read_header(opts.infp, opts.sep)
    cnames = cleangwas.parse_header(header, opts, logger, cleangwas.default_cnames)
    insnps = pd.read_csv(opts.include, header=None, names=["RS"], sep='\s+')
    insnps['RS'] = insnps['RS'].str.upper()
    total_df = pd.DataFrame()
    converter = cleangwas.get_converter(cnames)
    for chunk in pd.read_csv(opts.infp, sep=opts.sep, header=0, names=cnames, dtype=str, iterator=True,
                             chunksize=2000000):
        for k,v in converter.items(): chunk[k] = chunk[k].apply(v)
        if 'RS' in chunk.columns: 
            chunk['RS'] = chunk['RS'].str.upper()
            chunk = chunk[chunk['RS'].isin(insnps['RS'])]
            if not chunk.empty:
                total_df = total_df.append(chunk, sort=False)
        else:
            logger.logger.warning('Can not identify RS column for {0}\nExiting'.format(opts.infp))
            exit(1)
    if not total_df.empty:
        content = total_df.to_csv(path_or_buf=None, sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
        opts.infp = StringIO(content)
        df = cleangwas.qc(opts, cnames, logger)
        df = cleangwas.selectSNP(opts, df, logger)
        df = cleangwas.resort_col(df)
        df = cleangwas.truncate(opts, df)
    else:
        logger.logger.info('No common SNPs among exposure and outcome.')
        df = pd.DataFrame()
    try:
        os.remove('include')
    except :
        pass
    return df


def intersect(opts, dfs:tuple, logger):
    col_names = opts.identifier.strip().split(":")
    ids = list()
    for x in dfs:
        a1_a2 = set(x.apply(func=lambda series: ':'.join([series[x] for x in col_names]), axis=1))
        col_names = opts.identifier.strip().replace('1', '3').replace('2', '1').replace('3', '2').split(":")
        a2_a1 = set(x.apply(func=lambda series: ':'.join([series[x] for x in col_names]), axis=1))
        a1_a2_reverse = {ele.replace('A', '[').replace('T', 'A').replace('[', 'T').replace('C','[').replace('G','C').replace('[', 'G') for ele in a1_a2}
        a2_a1_reverse = {ele.replace('A', '[').replace('T', 'A').replace('[', 'T').replace('C', '[').replace('G', 'C').replace('[', 'G') for ele in a2_a1}
        ids.append(a1_a2 | a2_a1 | a1_a2_reverse | a2_a1_reverse)
    ids = reduce(lambda x,y : x & y, ids)
    ids = {x.replace('A', '').replace('T', '').replace('C', '').replace('G','').strip().strip(":") for x in ids}
    logger.logger.info('{0} common SNPs across {1} files.'.format(len(ids), len(dfs)))
    if len(ids) == 0:
        return [pd.DataFrame() for x in range(len(dfs))]
    col_names = opts.identifier.strip().replace('A1:A2','').strip(':').split(":")
    newdfs = list()
    for x in dfs:
        x['id____'] = x.apply(func=lambda series: ':'.join([series[x] for x in col_names]), axis=1)
        x = x[x['id____'].isin(ids)]
        x = x.drop(columns=['id____'])
        newdfs.append(x)
    return newdfs

def mr_test(expdf, outcdf, opts, logger):
    logger.logger.info('>>>>>>>>>>MR analyses for {0} and {1} start<<<<<<<<<<'.format(opts.exp_name, opts.out_name))
    expdf['RS'] = expdf['RS'].str.lower()
    outcdf['RS'] = outcdf['RS'].str.lower()
    expdf.to_csv('exposure.txt', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    outcdf.to_csv('outcome.txt', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    pleiotropy(opts,logger)
    mr_result = pd.DataFrame({'outcome':[],'exposure':[],'method':[],'nsnp':[],'b':[],'se':[], 'pval':[]})
    mr_heter = pd.DataFrame();mr_pleio = pd.DataFrame();mr_data = pd.DataFrame();mr_pleiosnp=pd.DataFrame()
    r_errf = tempfile.TemporaryFile(mode='w+', encoding='utf-8')
    if opts.heterogeneity != 'none':
        if os.path.exists('mr.pleiosnp'):
            mr_pleiosnp = pd.read_csv('mr.pleiosnp', sep='\s+', header=0, dtype=str)
            try:
                os.remove('mr.pleiosnp')
            except:
                pass
    if 'mr_gsmr' in opts.mr_method:
        if expdf.shape[0] < 5:
            logger.logger.info('No enough SNPs for GSMR analysis.')
        else:
            logger.logger.info('Notes: make sure that SNP frequency are accurate for GSMR analysis')
            refdir = os.path.split(os.path.realpath(__file__))[0]
            sysstr = system_platform()
            gctascript = os.path.join(refdir, 'refsource', 'gcta_{0}'.format(sysstr), 'gcta64')
            with open('gsmr_ref_data.txt', 'w') as writer:
                writer.write('\n'.join([os.path.join(refdir, 'refsource', '1000G_EUR_Phase3_plink', '1000G.EUR.QC.{0}'.format(str(i))) for i in range(1,23)]))
            with open('gsmr_exposure.txt', 'w') as writer:
                writer.write('\t'.join([opts.exp_name, 'exposure.txt'])+'\n')
            with open('gsmr_outcome.txt', 'w') as writer:
                writer.write('\t'.join([opts.out_name, 'outcome.txt'])+'\n')
            try:
                commands = '{0} --mbfile gsmr_ref_data.txt --gsmr-file gsmr_exposure.txt gsmr_outcome.txt --gsmr-direction 0 --out gsmr_result_pleio --heidi-thresh 0  --effect-plot --diff-freq 1 --clump-r2 1 --gwas-thresh 1 --gsmr-snp-min 5'.format(gctascript)
                content1 = subprocess.check_output(commands,
                                             stderr=None, shell=True, bufsize=0)
                '''gsmr_result = pd.read_csv('gsmr_result_pleio.gsmr', sep='\s+', header=0,dtype=str,
                                          names=['exposure', 'outcome', 'b', 'se', 'pval', 'nsnp', 'global_heidi_outlier'])'''
                gsmr_result = pd.read_csv('gsmr_result_pleio.gsmr', sep='\s+', header=0,dtype=str,
                                          names=['exposure', 'outcome', 'b', 'se', 'pval', 'nsnp'])
                gsmr_result['method'] = 'gsmr' * gsmr_result.shape[0]
                gsmr_result = gsmr_result[['outcome', 'exposure', 'method', 'nsnp', 'b', 'se', 'pval']].copy()

                mr_result = pd.concat([mr_result, gsmr_result], axis=0)
                #logger.logger.info(content1.decode('utf-8'))
                if opts.heterogeneity == 'none':
                    r_errf.seek(0)
                    r_errf.truncate()
                    content2 = subprocess.check_output('{0} --mbfile gsmr_ref_data.txt --gsmr-file gsmr_exposure.txt gsmr_outcome.txt --gsmr-direction 0 --out gsmr_result_nonpleio --heidi-thresh 0.01  --effect-plot --diff-freq 1 --clump-r2 1 --gwas-thresh 1 --gsmr-snp-min 5'.format(gctascript),
                                                 stderr=r_errf, shell=True, bufsize=0)
                    #logger.logger.info(content2.decode('utf-8'))
                    gsmr_result = pd.read_csv('gsmr_result_nonpleio.gsmr', sep='\s+', header=0,
                                              names=['exposure', 'outcome', 'b', 'se', 'pval', 'nsnp'], dtype=str)
                    gsmr_result['method'] = 'gsmr-outlier-correction' * gsmr_result.shape[0]
                    gsmr_result = gsmr_result[['outcome', 'exposure', 'method', 'nsnp', 'b', 'se', 'pval']].copy()
                    mr_result = pd.concat([mr_result, gsmr_result], axis=0)
            except subprocess.CalledProcessError as e:
                r_errf.flush()
                r_errf.seek(0)
                logger.logger.info('Error info:\n{0}'.format(r_errf.read()))
                logger.logger.info("GSMR analysis Failed\nSkiping...")
                r_errf.seek(0)
                r_errf.truncate()
        try:
            os.remove('gsmr_ref_data.txt')
            os.remove('gsmr_exposure.txt')
            os.remove('gsmr_outcome.txt')
            os.remove('gsmr_result_nonpleio.gsmr')
            os.remove('gsmr_result_pleio.gsmr')
            os.remove('gsmr_result_nonpleio.log')
            os.remove('gsmr_result_pleio.log')
            os.remove('gsmr_result_nonpleio.pleiosnps')
        except:
            pass
    if set(["mr_wald_ratio", "mr_two_sample_ml mr_egger_regression", "mr_egger_regression_bootstrap", "mr_simple_median", "mr_weighted_median",
        "mr_penalised_weighted_median", "Penalised weighted median", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre" ,"mr_ivw_fe", "mr_simple_mode",
        "mr_weighted_mode", "mr_weighted_mode_nome", "mr_simple_mode_nome" ,"mr_raps", "mr_sign", "mr_uwr"]) & set(opts.mr_method):
        try:
            #写提示，SNP少某些MR方法无法实现
            output = subprocess.check_output('Rscript mr_test.R {0} {1} {2} {3} {4} {5}'.format('exposure.txt', 'outcome.txt', opts.exp_name, opts.out_name, ','.join(set(opts.mr_method)), 'no' if opts.heterogeneity in ['mr_presso','mr_gsmr'] else 'yes'),
                                             stderr=r_errf, shell=True, bufsize=0)
            #logger.logger.info(output.decode('utf-8'))
            if os.path.exists('mr.result'):
                twosamplemr_result = pd.read_csv('mr.result', sep='\s+', header=0,dtype=str)
                mr_result = pd.concat([mr_result, twosamplemr_result], axis=0)
                os.remove('mr.result')
            if os.path.exists('mr.heter'):
                twosamplemr_heter = pd.read_csv('mr.heter', sep='\s+', header=0,dtype=str)
                mr_heter = pd.concat([mr_heter, twosamplemr_heter], axis=0)
                os.remove('mr.heter')
            if os.path.exists('mr.pleio'):
                twosamplemr_pleio = pd.read_csv('mr.pleio', sep='\s+', header=0,dtype=str)
                mr_pleio = pd.concat([mr_pleio,twosamplemr_pleio], axis=0)
                os.remove('mr.pleio')
            if os.path.exists('mr.data'):
                twosamplemr_data = pd.read_csv('mr.data',  sep='\s+', header=0,dtype=str)
                mr_data = pd.concat([mr_data, twosamplemr_data], axis=0)
                os.remove('mr.data')
        except subprocess.CalledProcessError as e:
            r_errf.flush()
            r_errf.seek(0)
            logger.logger.info('Error info:\n{0}'.format(r_errf.read()))
            logger.logger.info("MR analysis Failed.\nSkiping...")
        except FileNotFoundError as e1:
            pass

        try:
            os.remove('exposure.txt')
            os.remove('outcome.txt')
        except:
            pass
    r_errf.close()
    if mr_result.empty:
        logger.logger.info('\n=====No MR results=====\n')
    else:

        logger.logger.info('\n=====Main MR results=====\n')
        logger.logger.info(str(mr_result))
    if mr_heter.empty:
        logger.logger.info('\n=====No heterogeneity test results=====\n')
    else:
        logger.logger.info('\n=====Main heterogeneity test for MR results=====\n')
        logger.logger.info(str(mr_heter))
    if mr_pleio.empty:
        logger.logger.info('\n=====No Egger pleiotropy test result=====\n\n')
    else:
        logger.logger.info('\n--Egger pleiotropy: Intercetpt--\n\n')
        logger.logger.info(str(mr_pleio))
    logger.logger.info('\n\n>>>>>>>>>>MR analysis for {0} and {1} finished!<<<<<<<<<<\n'.format(opts.exp_name, opts.out_name))
    return mr_result,mr_heter,mr_pleio,mr_data,mr_pleiosnp



def mr4general(opts, logger):
    abspath = os.path.split(os.path.realpath(__file__))[0]
    expdf = read_exp(opts, logger)
    clumping(opts, expdf, abspath, logger)
    outcdf = read_outc(opts, logger)
    expdf, outcdf = intersect(opts, tuple([expdf, outcdf]), logger)
    if outcdf.empty:
        logger.logger.info("As no common SNPs between {0} and {1}, MR analysis for this trait pair terminated!".format(opts.exp_name, opts.out_name))
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp = mr_test(expdf, outcdf,  opts, logger)
    return mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp

def mr4gm(opts, logger):
    opts = copy.deepcopy(opts)
    opts.infp = opts.exposure_gm
    header = cleangwas.read_header(opts.infp, opts.sep)
    cnames = cleangwas.parse_header(header, opts, logger, cleangwas.default_cnames)
    all_df = pd.read_csv(opts.exposure_gm, header=0, names=cnames, sep='\s+', dtype=str)
    all_df['RS'].to_csv('include', header=False, sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    out_df = read_outc(opts, logger)
    mr_results = []; mr_heters=[]; mr_pleios=[]; mr_datas=[];mr_pleiosnps = []
    for phe in [x for x in all_df.columns if x in [ 'Phenotype']]:
    #for phe in [x for x in all_df.columns if x in [ 'Phenotype', 'Category']]:
        for x in set(all_df.loc[:,phe]):
            all_df[all_df[phe] == x].copy().to_csv(path_or_buf='.exposure_gm.txt', sep='\t', na_rep='NA', float_format='%g',
                                                   encoding='utf-8', index=False)
            out_df.to_csv(path_or_buf='.outcome_gm.txt', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
            opts.exposure = '.exposure_gm.txt'
            opts.outcome = '.outcome_gm.txt'
            opts.exp_name = x
            opts.exposure_gm = None
            mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp = mr4general(opts, logger)
            mr_results.append(mr_result); mr_heters.append(mr_heter); mr_pleios.append(mr_pleio); mr_datas.append(mr_data);mr_pleiosnps.append(mr_pleiosnp)
            try:
                os.remove('.exposure_gm.txt')
                os.remove('.outcome_gm.txt')
            except :
                pass
    if len([x for x in mr_results if not x.empty]) != 0:

        mr_results = pd.concat([x for x in mr_results if not x.empty], axis=0, ignore_index=True,
              sort=False, copy=True)
    if len([x for x in mr_heters if not x.empty]) != 0:
        mr_heters = pd.concat([x for x in mr_heters if not x.empty], axis=0, ignore_index=True,
              sort=False, copy=True)
    if len([x for x in mr_pleios if not x.empty]) !=0:
        mr_pleios = pd.concat([x for x in mr_pleios if not x.empty], axis=0, ignore_index=True,
                          sort=False, copy=True)
    if len([x for x in mr_datas if not x.empty]) != 0:
        mr_datas = pd.concat([x for x in mr_datas if not x.empty], axis=0, ignore_index=True,
              sort=False, copy=True)
    if opts.heterogeneity!='none' and len([x for x in mr_pleiosnps if not x.empty]) !=0:
        mr_pleiosnps = pd.concat([x for x in mr_pleiosnps if not x.empty], axis=0, ignore_index=True,
                             sort=False, copy=True)
    else:
        mr_pleiosnps = pd.DataFrame()
    return mr_results,mr_heters, mr_pleios, mr_datas, mr_pleiosnps
def write2file(mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp, opts, logger):
    mr_result.to_csv(opts.outfp+'.mr_result', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    mr_heter.to_csv(opts.outfp+'.mr_heter', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    mr_pleio.to_csv(opts.outfp+'.mr_pleio', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    mr_data.to_csv(opts.outfp+'.mr_data', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    if opts.heterogeneity != 'none':
        mr_pleiosnp.to_csv(opts.outfp+'.mr_pleiosnp', sep='\t', na_rep='NA', float_format='%g', encoding='utf-8', index=False)
    logger.logger.info('''Results have been written into file {0}.mr_result, {0}.mr_heter, {0}.mr_pleio, {0}.mr_data, {0}.mr_pleiosnp.'''.format(opts.outfp))


logger = Logger()

@timeit(logger)
def init(opts, logger):
    if opts.exposure:
        mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp = mr4general(opts,logger)
    else:
        mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp = mr4gm(opts,logger)
    write2file(mr_result, mr_heter, mr_pleio, mr_data, mr_pleiosnp, opts, logger)
if __name__ == '__main__':
    opts = parse_arg()
    init(opts, logger)
    logger.logger.info('End!')
    logger.rename(opts.outfp)
    sys.exit(0)


