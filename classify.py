import gzip
import sys
import logging
import math

def read_inv_list(inv_list):
    invinfo=open(inv_list, 'r').read().split('\n')[:-1]
    if invinfo[0][0]=='#':
        invinfo=invinfo[1:]
    inv_list={}
    for inv in invinfo:
        inv = inv.split('\t')
        inv_list[inv[0]] = inv[1]
    return inv_list


def distance(a,b):
    if len(a)!=len(b):
        print('length of two position do not match.')
        print(len(a),len(b))
        quit()
    sumdis=0
    for i in range(len(a)):
        sumdis+=(float(a[i])-b[i])**2

    sumdis=sumdis**0.5
    return sumdis


def genotype(inv, profile_path, vcffile, sample_idx):
    # read tag snp info from profile
    tagsnp=open(profile_path,'r').read().split('\n')[:-1]
    chrom=tagsnp[0].split('\t')[1]
    constant=[float(c.split('\t')[5])*abs(float(c.split('\t')[5])) for c in tagsnp]
    allpos=[int(c.split('\t')[2]) for c in tagsnp]
    min_tag_pos=min(allpos)
    max_tag_pos=max(allpos)

    print('Start genotyping for:'+inv)
    samplesnp={}
    sampleid=[]
    snp_genotype={}
    numsamples=len(sample_idx)-9
    for sample in sample_idx[9:]:
        snp_genotype[sample]={}
        snp_genotype[sample]['Homo']=0
        snp_genotype[sample]['Hetero']=0

    sampleid=sample_idx[9:]
    for i in range(len(sampleid)):
        samplesnp[sampleid[i]] = [-1.0*abs(c) for c in constant]

    # get list of SNP calls from VCF
    try:
        allsnpcall = vcffile.query(chrom, min_tag_pos-1000, max_tag_pos+1000)
    except:
        print('Warning: Failed to extract SNP for '+inv+', skipping...')
        return []

    detected_tagsnp=0
    for snpcall in allsnpcall:
        # skip non-SNPs
        if len(snpcall[3])!=1 or len(snpcall[4])!=1:
            continue

        gtindex=snpcall[8].split(':').index("GT")
        istagsnp=0
        for tagidx in range(len(tagsnp)):
            snp=tagsnp[tagidx].split('\t')
            if snpcall[1]==snp[2] and snp[3]==snpcall[3] and snp[4]==snpcall[4]:
                istagsnp=1
                break

        # skip non-tag SNPs
        if istagsnp ==0:
            continue
        detected_tagsnp+=1

        # find GT position for each sample
        gtinfo=snpcall[9:]
        gtinfo=[c.split(':')[gtindex] for c in gtinfo]

        for i in range(numsamples):
            if '1' in gtinfo[i] or '2' in gtinfo[i]:
                samplesnp[sampleid[i]][tagidx] = abs(samplesnp[sampleid[i]][tagidx])
            if gtinfo[i] in ['.|.','./.']:
                samplesnp[sampleid[i]][tagidx] = 0
            if gtinfo[i] in ['1/1','1|1']:
                snp_genotype[sampleid[i]]['Homo']+=1
            elif gtinfo[i] in ['1/0','0/1','1|0','0|1']: 
                snp_genotype[sampleid[i]]['Hetero']+=1
    f=open('matrics_inv_tag-'+inv,'w')
    for samplei in range(len(samplesnp)):
        for tagi in range(len(samplesnp[sampleid[samplei]])):
            f.write(str(samplesnp[sampleid[samplei]][tagi])+'\t')
        f.write('\n')
    f.close()
    print('Found ',detected_tagsnp,'/',len(tagsnp), ' tag SNPs for ',inv)

    # calculate distance to pos and neg center
    inv_calls=[]
    positivepos=constant
    for samplei in range(len(samplesnp)):
        dipositive=distance(samplesnp[sampleid[samplei]],positivepos)
        dinegative=distance(samplesnp[sampleid[samplei]], [-1*c for c in positivepos])

        if snp_genotype[sampleid[samplei]]['Hetero']>snp_genotype[sampleid[samplei]]['Homo']:
            invgenotype='Hetero'
        else:
            invgenotype='Homo'
        inv_calls+=[[sampleid[samplei],dipositive,dinegative,invgenotype,snp_genotype[sampleid[samplei]]['Homo'],snp_genotype[sampleid[samplei]]['Hetero']]]

    return inv_calls


def write_output(inv_results, output_prefix, min_score):
    outfile=open(output_prefix,'w')
    for inv in inv_results:
        for invcall in inv_results[inv]:
            score=-10*math.log10((invcall[1]+0.001)/(invcall[2]+0.001))
            if score>=min_score:
                score=str(round(score*10)/10)
                outfile.write(invcall[0]+'\t'+inv+'\t'+invcall[3]+'\t'+score+'\t'+str(invcall[1])+'\t'+str(invcall[2])+'\t'+str(invcall[4])+'\t'+str(invcall[5])+'\n')
    outfile.close()

