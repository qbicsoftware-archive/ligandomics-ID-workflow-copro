from CTDopts.CTDopts import _InFile, CTDModel, args_from_file
import sys
import os
import subprocess
import re
import pandas as pd
from itertools import combinations


pattern = re.compile('Q\w{4}[0-9]{3}[a-zA-Z]\w')

wf_dir = sys.argv[1]
ctd_params = args_from_file(wf_dir + '/WORKFLOW-CTD')
ctd_files = args_from_file(wf_dir + '/IN-FILESTOSTAGE')

data_path = '%s/data/' % wf_dir
result_path = '%s/result/' % wf_dir
log_path = '%s/logs/' % wf_dir
db_path = os.path.join(wf_dir, 'ref')

mzmlFiles = []

# mzml files
for filePath in ctd_files['Mass Spectrometry Data']:
    fileName = filePath.split('/')[-1]
    mzmlFiles.append('%s%s' % (data_path, fileName))

# Parameters
fmt = float(ctd_params['fmt'])
pmt = float(ctd_params['pmt'])
fbo = float(ctd_params['fbo'])
fdr = float(ctd_params['fdr'])
num_hits = int(ctd_params['noh'])
dmr = ctd_params['dmr']

logfilename = 'ligandomicsID_v2_0_coPro_workflow.logs'
logfile = open(logfilename, 'w')

if ctd_files['db'] != '':
    fasta_path = os.path.join(db_path, ctd_files['db'].split('/')[-1])
else:
    fasta_path = os.path.join(data_path, ctd_files['Individualized Reference'].split('/')[-1])

fasta_decoy_path = os.path.join(data_path, 'proteinswithdecoys.fasta')

commandDecoy = 'DecoyDatabase  -in {i} -out {o} -decoy_string XXX -decoy_string_position prefix'.format(i=fasta_path, o=fasta_decoy_path)
subprocess.call(commandDecoy.split(), stderr=logfile, stdout=logfile)

idFiles = []


#mzmlFiles = ['/lustre_cfc/qbic/chris/playground/ligandomicsID/data/spectra_comet.mzML', '/lustre_cfc/qbic/chris/playground/ligandomicsID/data/spectra_comet_2.mzML', '/lustre_cfc/qbic/chris/playground/ligandomicsID/data/spectra_comet_3.mzML']
#idFiles = ['/lustre_cfc/qbic/chris/playground/ligandomicsID/data/spectra_comet.idXML', '/lustre_cfc/qbic/chris/playground/ligandomicsID/data/spectra_comet_2.idXML', '/lustre_cfc/qbic/chris/playground/ligandomicsID/data/spectra_comet_3.idXML']

mzmlFiles = ['1.mzML', '2.mzML', '3.mzML']

for mzml in mzmlFiles:
    if mzml.endswith('.gz'):
        logfile.write("Extracting gzipped content... \n")
        cmd = "gzip -d {f}".format(f=mzml)
        os.system(cmd)
        mzml = mzml.replace('.gz', '')

    idPath = mzml.replace('mzML', 'idXML')

    identifier = mzml.split('/')[-1].split('.')[0]

    if ctd_params['centroided'] == 'false':
        pickpeakcommand = 'PeakPickerHiRes -in {i} -out {o} -threads 20 -algorithm:ms_levels 1'
        subprocess.call(pickpeakcommand.format(i=mzml, o=mzml).split(),stderr=logfile, stdout=logfile)

    ##### NEW  #####

    commandComet = 'CometAdapter -in {i} -out {o} -threads 20 -database {d} -precursor_mass_tolerance {pmt} -fragment_bin_tolerance {fmt} -fragment_bin_offset {fbo} -num_hits {n} -digest_mass_range {dmr}'.format(i=mzml, o=idPath, d=fasta_decoy_path, pmt=pmt, fmt=fmt, fbo=fbo, n=num_hits, dmr=dmr) 
    subprocess.call(commandComet.split() + ["-fixed_modifications", "Carbamidomethyl (C)", "-variable_modifications", "Oxidation (M)", "-enzyme", "unspecific cleavage"],stderr=logfile, stdout=logfile)

    peptideIndexer = 'PeptideIndexer -in {f} -out {o} -threads 20 -fasta {d} -decoy_string XXX -enzyme:specificity none -enzyme:name '.format(f=idPath, o=idPath, d=fasta_decoy_path)
    subprocess.call(peptideIndexer.split()  + ["unspecific cleavage"],stderr=logfile, stdout=logfile)

    idFiles.append(idPath)



############# insert id based map transformation

alignedMaps=[a.replace('.idXML', '_aligned.idXML') for a in idFiles]
trafoMaps=[a.replace('.idXML', '.trafoXML') for a in idFiles]
MapAlignerID = 'MapAlignerIdentification -in {f} -out {o} -trafo_out {t} -threads 20'.format(f=" ".join(idFiles), o=" ".join(alignedMaps), t=" ".join(trafoMaps))
subprocess.call(MapAlignerID.split(), stderr=logfile, stdout=logfile)

for Map in trafoMaps:
    mzml = mzmlFiles[trafoMaps.index(Map)]
    trafoMZML=mzml.replace('.mzML', '_trafo.mzML')
    MapTransformer = 'MapRTTransformer -in {f} -out {o} -trafo_in {m} -threads 20'.format(m=Map, f=mzml, o=trafoMZML)
    subprocess.call(MapTransformer.split(), stderr=logfile, stdout=logfile)

idFiles=alignedMaps
#############



### predict hits of fitting length and calc FDR and PEP
#FDR calc

## FileMerging
#### change f and o

identifer_for_files = idFiles[0].split('/')[-1]

idresult = os.path.join(result_path, identifer_for_files.replace('.idXML','_merged.idXML'))
FileMerger = 'IDMerger -in {f} -out {o} -threads 20 -annotate_file_origin'.format(f=' '.join(idFiles),o=idresult)
subprocess.call(FileMerger.split(),stderr=logfile, stdout=logfile)

idresult_fdr = os.path.join(result_path, identifer_for_files.replace('.idXML', '_merged_fdr.idXML'))
falseDiscovery = 'FalseDiscoveryRate -in {f} -out {o} -threads 20 -algorithm:add_decoy_peptides -algorithm:use_all_hits'.format(f=idresult,o=idresult_fdr)
subprocess.call(falseDiscovery.split(),stderr=logfile, stdout=logfile)

### extract Percolator Features with PSMFeatureExtractor on Mathias TopPerc branch
idresult_fdr_psm = os.path.join(result_path, identifer_for_files.replace('.idXML', '_merged_fdr_psm.idXML'))
PSMFeat = 'PSMFeatureExtractor -in {f} -out {o} -threads 20'.format(f=idresult_fdr,o=idresult_fdr_psm)
subprocess.call(PSMFeat.split(),stderr=logfile, stdout=logfile)

### run Percolator with PercolatorAdapter on Mathias TopPerc branch
idresult_perc = os.path.join(result_path, '{}_merged_perc.idXML'.format(identifer_for_files))
Percolator = 'PercolatorAdapter -in {f} -out {o} -decoy-pattern XXX -debug 10 -threads 20 -enzyme no_enzyme -trainFDR 0.05 -testFDR 0.05'.format(f=idresult_fdr_psm,o=idresult_perc)
subprocess.call(Percolator.split(),stderr=logfile, stdout=logfile)

#filter by provided FDR value
idresult_filtered = os.path.join(result_path, '{}_merged_perc_fdr_filtered.idXML'.format(identifer_for_files))
idFilter = 'IDFilter  -in {f} -out {o} -score:pep {m} -remove_decoys -threads 20'.format(f=idresult_fdr_psm, o=idresult_filtered, m=fdr)
subprocess.call(idFilter.split(),stderr=logfile, stdout=logfile)

#IDRipper
#idresult_folder = os.path.join(result_path, '{}'.format(identifer_for_files))
#IDRipper = 'IDRipper  -in {f} -out {o} -threads 20'.format(f=idresult_filtered, o=idresult_folder)
#subprocess.call(IDRipper.split(),stderr=logfile, stdout=logfile)

#IDMerger
features = []
files = [i for i in idFiles]

for file in files:
   internal=idresult_filtered
   ### map percolator refined and FDR filtered ids onto features
   mergeresult = os.path.join(result_path, file.replace('idXML', 'featureXML').split('/')[-1])
   features.append(mergeresult)
   IDMapper = 'FeatureFinderIdentification -in {f} -id {i} -threads 20 -out {o}'.format(f=file.replace('_aligned.idXML','_trafo.mzML'),i=internal,o=mergeresult)
   subprocess.call(IDMapper.split(),stderr=logfile, stdout=logfile)

## FeatureLinking
mergeresult = os.path.join(result_path, idFiles[0].split('/')[-1].replace('idXML', 'consensusXML'))
FeatLinker = 'FeatureLinkerUnlabeledKD -in {f} -out {o} -threads 20'.format(f=' '.join(features), o=mergeresult)
subprocess.call(FeatLinker.split(),stderr=logfile, stdout=logfile)

## IDConflictResolver
ConfSolver = 'IDConflictResolver -in {f} -out {o}'.format(f=mergeresult,o=mergeresult)
subprocess.call(ConfSolver.split(),stderr=logfile, stdout=logfile)


########################################################################################

convert = 'TextExporter -in {f} -out {o} -id:add_hit_metavalues 0'.format(f=mergeresult, o=mergeresult.replace('.consensusXML', '.csv'))
subprocess.call(convert.split(),stderr=logfile, stdout=logfile)

########################################################################################

op=open(mergeresult.replace('.consensusXML', '.csv'))
opr=op.readlines()
op.close()

df = []

for i, r in enumerate(opr):
    if r.startswith('#PEPTIDE'):
        header = r.strip().split('\t')[1:] + opr[i - 1].strip().split('\t')[1:]
    if r.startswith('PEPTIDE'):
        if not opr[i-1].startswith('PEPTIDE'):
            df.append(r.strip().split('\t')[1:] + opr[i - 1].strip().split('\t')[1:])

df=pd.DataFrame(df)
df.columns=header
df.to_csv(mergeresult.replace('.consensusXML', '_edit.csv'))

seqs=df['sequence'].values.tolist()
seqs_new=[]
for seq in seqs:
    m = re.findall("\(\w+\)", seq)
    for i in m:
        seq = seq.replace(i, "")
    seqs_new.append(seq)

seqs_new=list(set(seqs_new))
df['sequence']=df['sequence'].str.replace("\(\w+\)","")

Final_df={}
Final_df_header=['fdr', 'xcorr', 'deltacn', 'median_intensity']
for seq in seqs_new:
    row=df[df['sequence']==seq].sort('score').iloc[[0]]
    Final_df[seq]=[row['score'].values.tolist()[0],row['MS:1002252'].values.tolist()[0],row['MS:1002253'].values.tolist()[0],row['intensity_cf'].values.tolist()[0]]
Final_df=pd.DataFrame(Final_df).transpose()
Final_df.columns=Final_df_header
Final_df.to_csv(mergeresult.replace('.consensusXML', '_extracted.csv'))

logfile.close()
subprocess.call(["mv", logfilename, log_path])
