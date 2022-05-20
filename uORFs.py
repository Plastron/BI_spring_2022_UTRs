import urllib

#complementary table
compl = {"A":"T", "T":"A", "G":"C", "C":"G"}

#genetic code table
gen_code = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}



probs = {}

probs_file = open("D:/bioinformatic_institution/uORF_project/probs.txt", 'r')

#converting trinucleotide mutation probabilities from the file into the dict {NXN: [Y, probability]}
for each_line in probs_file:
    inner_line = each_line.split(" ")
    if inner_line[0] in probs:
        probs[inner_line[0]][inner_line[1]] = float(inner_line[2][:-1])
    else:
        probs[inner_line[0]] = {inner_line[1]: float(inner_line[2][:-1])}

probs_file.close()

#the uORFs file
bed_file = open("D:/bioinformatic_institution/uORF_project/Our_AD+AR_uORFs_v3.bed", 'r')

#file with resulting probabilities
res_probs = open("D:/bioinformatic_institution/uORF_project/res_probs_2.txt", 'w')

counter = 0
for line in bed_file:
    if (counter != 0) and ("StartUpstreamMainFrameStart|StopUpstreamMainFrameStop|StopUpstreamMainFrameStart" in line): #the second condition is to choose uORFs which do not overlap with the main ORF
        inner_line = line.split()
        sequence = urllib.request.urlopen("http://togows.org/api/ucsc/hg38/" + inner_line[0] + ":" + str(int(inner_line[1]) + 1) + "-" + inner_line[2] + ".fasta").read() # the URL request to obtain nucleotide sequene
        sequence = sequence.decode()
        sequence = sequence.split('\n') # the sequence is in FASTA format with several lines. we need to cancatenate it in a single string
        name = sequence[0]
        whole_seq = "".join(sequence[1:])

        # now we need to obtain the coding sequence, so we splice the introns
        exons_only = ""
        starts = [int(i) for i in inner_line[11].split(',')[:-1]] #the coordinates of exon starts
        length = [int(i) for i in inner_line[10].split(',')[:-1]] #the lengths of exons
        if len(starts) > 1:
            for i in range(len(starts)):
                exons_only += whole_seq[starts[i]:starts[i]+length[i]]
        else:
            exons_only = whole_seq

        if inner_line[5] == "-": #if the uORF is on the 'minus' strand, we need to obtain its complementary sequence for analysis
            exons_only = "".join([compl[i] for i in exons_only[::-1]])

        #next we calculate the probabilities of mutations
        syn_prob = 0
        mis_prob = 0
        non_prob = 0
        if len(exons_only) <= 6: #there are uORFs of length 6, which means only the start and stop codons are in it. We don't use these for our analysis
            syn_prob = "NA"
            mis_prob = "NA"
            non_prob = "NA"
        else:
            for i in range(3,len(exons_only)-4,3): #scan through each codon in uORF, not including start and stop codons
                codon = exons_only[i:i+3]
                acid = gen_code[codon] #aminoacid encoded by the codon

                #next is the implementation of the method illustrated above
                #we separately look through the first, the second and the sird positions of the cdon
                #for each substitution in the position we calculate the outcome (synonymous/missense/nonsense) and obtain probability
                #then we add the probability to the outcome for the uORF
                first_tri = exons_only[i-1:i+2]
                try:
                    first_probs = probs[first_tri]
                #our probability table contains only the halve of all possible trinucleotides, assuming that, for example, AAA and TTT are the same in terms of probabilities
                except:
                    inner_codon = "".join([compl[i] for i in first_tri[::-1]])
                    inner_probs = probs[inner_codon]
                    first_probs = {compl[i]:inner_probs[i] for i in inner_probs}
                #substitute each nucleotide in the first position of the triplet to the other three, look for the outcome and add the probability to the total probability for the outcome
                for each_sub in first_probs:
                    res_tri = each_sub + codon[1:]
                    if gen_code[res_tri] == acid:
                        syn_prob += first_probs[each_sub]
                    elif gen_code[res_tri] == "*":
                        non_prob += first_probs[each_sub]
                    else:
                        mis_prob += first_probs[each_sub]

                #the same operation but for the second position
                second_tri = codon
                try:
                    second_probs = probs[second_tri]
                except:
                    inner_codon = "".join([compl[i] for i in second_tri[::-1]])
                    inner_probs = probs[inner_codon]
                    second_probs = {compl[i]:inner_probs[i] for i in inner_probs}
                for each_sub in second_probs:
                    res_tri = codon[0] + each_sub + codon[2]
                    if gen_code[res_tri] == acid:
                        syn_prob += second_probs[each_sub]
                    elif gen_code[res_tri] == "*":
                        non_prob += second_probs[each_sub]
                    else:
                        mis_prob += second_probs[each_sub]

                #and for the sird position
                sird_tri = exons_only[i+1:i+4]
                try:
                    sird_probs = probs[sird_tri]
                except:
                    inner_codon = "".join([compl[i] for i in sird_tri[::-1]])
                    inner_probs = probs[inner_codon]
                    sird_probs = {compl[i]:inner_probs[i] for i in inner_probs}
                for each_sub in sird_probs:
                    res_tri = codon[:2] + each_sub
                    if gen_code[res_tri] == acid:
                        syn_prob += sird_probs[each_sub]
                    elif gen_code[res_tri] == "*":
                        non_prob += sird_probs[each_sub]
                    else:
                        mis_prob += sird_probs[each_sub]
        #write the chromosome, start, end, name, strand and final probabilities into the final file
        res_probs.write(inner_line[0] + " " + inner_line[1] + " " + inner_line[2] + " " + inner_line[3] + " " + inner_line[5] + " " + exons_only + " " + str(syn_prob)+ " " + str(mis_prob)+ " " + str(non_prob) + "\n")
    counter += 1
    print(counter)

bed_file.close()
res_probs.close()


import urllib.request as urllib2
from bs4 import *
from urllib.parse  import urljoin
url = "https://gnomad.broadinstitute.org/region/12-70243450-70243495?dataset=gnomad_r3"
r = urllib.request.urlopen(url)
print(r.full_url)

def crawl(pages, depth=None):
    indexed_url = [] # a list for the main and sub-HTML websites in the main website
    for i in range(depth):
        for page in pages:
            if page not in indexed_url:
                indexed_url.append(page)
                try:
                    c = urllib2.urlopen(page)
                except:
                    print( "Could not open %s" % page)
                    continue
                soup = BeautifulSoup(c.read())
                links = soup('a') #finding all the sub_links
                for link in links:
                    if 'href' in dict(link.attrs):
                        url = urljoin(page, link['href'])
                        if url.find("'") != -1:
                                continue
                        url = url.split('#')[0]
                        if url[0:4] == 'http':
                                indexed_url.append(url)
        pages = indexed_url
    return indexed_url


pagelist=["https://gnomad.broadinstitute.org/region/12-70243450-70243495?dataset=gnomad_r3"]
urls = crawl(pagelist, depth=1)
print(urls)

url = "https://gnomad.broadinstitute.org/region/12-70243450-70243495?dataset=gnomad_r3"
c = urllib2.urlopen(url)
soup = BeautifulSoup(c.read())
print(soup.prettify())








from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
import time


# configure webdriver
options = Options()
options.headless = True  # hide GUI
options.add_argument("--window-size=1920,1080")  # set window size to native GUI size
options.add_argument("start-maximized")  # ensure window is full-screen


# configure chrome browser to not load images and javascript
chrome_options = webdriver.ChromeOptions()
chrome_options.add_experimental_option("prefs", {"profile.managed_default_content_settings.images": 2})

driver = webdriver.Chrome(service=Service('D:/bioinformatic_institution/uORF_project/chromedriver_win32/chromedriver.exe'))

coords_file = open('D:/bioinformatic_institution/uORF_project/lost_uORFs.txt', 'r')

for line in coords_file:
    inner_line = line.split(" ")
    #form the URL request
    coords = inner_line[0][3:] + "-" + inner_line[1] + "-" + inner_line[2]
    url = "https://gnomad.broadinstitute.org/region/" + coords + "?dataset=gnomad_r3"

    driver.get(url)
    time.sleep(20) #it requires about 2 seconds for GnomAD to generate each page, so I was lazy and just gave it 20 seconds for a certain generation.

    for elm in driver.find_elements(By.TAG_NAME, "button"):
        if elm.text == 'Export variants to CSV':
            elm.click()

driver.close()


#analyse the gnomAD results
import os
dir_str = 'D:/bioinformatic_institution/uORF_project/gnomAD'
directory = os.fsencode(dir_str)

coords_dict = {}
#open each downloaded GnomAD file,
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    coord_str = filename.split("_v3.1.2_")[1].split("_2022_04_")[0]
    coords = coord_str.split('-')
    inner_gnomad_file = open(dir_str + '/' + filename, 'r')
    counter = 0
    for each_str in inner_gnomad_file:
        if counter > 0: #that's just to exclude the first string
            inner_var_str = each_str.split(',')
            if inner_var_str[1] != coords[2]: #that's to exclude the last position in gnomad file due to differences in gnomad and BED coordinates
                if coord_str in coords_dict:
                    coords_dict[coord_str].append(inner_var_str)
                else:
                    coords_dict[coord_str] = [inner_var_str]
        counter += 1
    inner_gnomad_file.close()


all_uORs_file = open('D:/bioinformatic_institution/uORF_project/Our_AD+AR_uORFs_v3.bed', 'r')

exons_coords = {}

#for each uORF write the exones coordinates to exclude intron variants in gnomAD file from calculations
counter = 0
for each_line in all_uORs_file:
    if counter > 0:
        inner_line = each_line.split()
        coords = inner_line[0][3:] + "-" + inner_line[1] + "-" + inner_line[2]
        inner_exons_coords = []
        starts = [int(i) for i in inner_line[11].split(',')[:-1]]
        length = [int(i) for i in inner_line[10].split(',')[:-1]]

        total_length = 0
        for i in range(len(starts)):
            inner_exons_coords.append([int(inner_line[1]) + starts[i], int(inner_line[1]) + starts[i] + length[i], total_length])
            total_length += length[i]

        exons_coords[coords] = inner_exons_coords
    counter += 1

all_uORs_file.close()
miss_counter = 0
correct_counter = 0

uORF_coords_file = open('D:/bioinformatic_institution/uORF_project/res_probs_2.txt', 'r')
uORF_coords_with_real_vars = open('D:/bioinformatic_institution/uORF_project/uORF_coords_with_real_vars.txt', 'w')

#now we need (1) for each of the gnomAd variant see if it's in exon
#(2) see the outcome of the variant for uORF product
for each_str in uORF_coords_file:
    real_syn = 0
    real_mis = 0
    real_non = 0
    inner_line = each_str.split(" ")
    coords = inner_line[0][3:] + "-" + inner_line[1] + "-" + inner_line[2]

    exons_intervals = exons_coords[coords]

    if coords in coords_dict and len(inner_line[5]) > 6:
        for each_var in coords_dict[coords]:
            if len(each_var[3]) == len(each_var[4]) == 1: #we only analyse SNP, no deletions r insertions as we didn't calculate the probabilities for them
                for each_exon_interval in exons_intervals:
                    var_position = int(each_var[1]) - 1
                    if (each_exon_interval[0] <= var_position) and (var_position < each_exon_interval[1]):
                        rel_var_pos = each_exon_interval[2] + var_position - each_exon_interval[0]
                        #if the uORF is on the "minus" strand, we invert the relative coordinates and look for complementary base
                        if inner_line[4] == "-":
                            rel_var_pos = len(inner_line[5]) - 1 - rel_var_pos
                            each_var[3] = compl[each_var[3]]
                            each_var[4] = compl[each_var[4]]

                        try:
                            #this block simply checks if the expected base from gnomAD and the base in our sequence are the same
                            if inner_line[5][rel_var_pos] != each_var[3]:
                                miss_counter += 1
                                print(str(rel_var_pos) + ' ' + inner_line[5] + ' ' + inner_line[4] + ' ' + each_var[3] + ' ' + inner_line[5][rel_var_pos])
                            else:
                                correct_counter += 1
                        except:
                            continue
                        try:
                            if rel_var_pos % 3 == 0:
                                ref_codon = inner_line[5][rel_var_pos:rel_var_pos + 3]
                                ref_amino = gen_code[ref_codon]
                                var_codon = each_var[4] + inner_line[5][rel_var_pos + 1:rel_var_pos + 3]
                                var_amino = gen_code[var_codon]
                                if var_amino == "*":
                                    real_non += 1
                                elif ref_amino == var_amino:
                                    real_syn += 1
                                else:
                                    real_mis += 1
                            elif (rel_var_pos - 1) % 3 == 0:
                                ref_codon = inner_line[5][rel_var_pos - 1:rel_var_pos + 2]
                                ref_amino = gen_code[ref_codon]
                                var_codon = inner_line[5][rel_var_pos - 1] + each_var[4] + inner_line[5][rel_var_pos + 1]
                                var_amino = gen_code[var_codon]
                                if var_amino == "*":
                                    real_non += 1
                                elif ref_amino == var_amino:
                                    real_syn += 1
                                else:
                                    real_mis += 1
                            else:
                                ref_codon = inner_line[5][rel_var_pos - 2:rel_var_pos + 1]
                                ref_amino = gen_code[ref_codon]
                                var_codon = inner_line[5][rel_var_pos -2:rel_var_pos] + each_var[4]
                                var_amino = gen_code[var_codon]
                                if var_amino == "*":
                                    real_non += 1
                                elif ref_amino == var_amino:
                                    real_syn += 1
                                else:
                                    real_mis += 1
                        except:
                            print(str(rel_var_pos) + ' ' + inner_line[5] + ' ' + inner_line[4])


        uORF_coords_with_real_vars.write(each_str[:-1] + " " + str(real_syn) + " " + str(real_mis) + " " + str(real_non) + "\n")


print(miss_counter) #the number of positions in which the GnomAD base differed from our sequence base.
print(correct_counter)#the number of positions in which the GnomAD base matched with our sequence base.
uORF_coords_file.close()
uORF_coords_with_real_vars.close()



sequence = urllib.request.urlopen("http://togows.org/api/ucsc/hg38/chr16:57245264-57245294.fasta").read()
sequence = sequence.decode()
sequence = sequence.split('\n')
name = sequence[0]
whole_seq = "".join(sequence[1:])
print(name)
print(whole_seq)


#reanalyse the gnomAD results: only SNP, + clinical significance + singletons

import os
dir_str = 'D:/bioinformatic_institution/uORF_project/gnomAD'
directory = os.fsencode(dir_str)

coords_dict = {}

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    coord_str = filename.split("_v3.1.2_")[1].split("_2022_04_")[0]
    coords = coord_str.split('-')
    inner_gnomad_file = open(dir_str + '/' + filename, 'r')
    counter = 0
    for each_str in inner_gnomad_file:
        if counter > 0:
            inner_var_str = each_str.split(',')
            if inner_var_str[1] != coords[2]:
                if coord_str in coords_dict:
                    coords_dict[coord_str].append(inner_var_str)
                else:
                    coords_dict[coord_str] = [inner_var_str]
        counter += 1
    inner_gnomad_file.close()


all_uORs_file = open('D:/bioinformatic_institution/uORF_project/Our_AD+AR_uORFs_v3.bed', 'r')

exons_coords = {}

counter = 0
for each_line in all_uORs_file:
    if counter > 0:
        inner_line = each_line.split()
        coords = inner_line[0][3:] + "-" + inner_line[1] + "-" + inner_line[2]
        inner_exons_coords = []
        starts = [int(i) for i in inner_line[11].split(',')[:-1]]
        length = [int(i) for i in inner_line[10].split(',')[:-1]]

        total_length = 0
        for i in range(len(starts)):
            inner_exons_coords.append([int(inner_line[1]) + starts[i], int(inner_line[1]) + starts[i] + length[i], total_length])
            total_length += length[i]

        exons_coords[coords] = inner_exons_coords
    counter += 1

all_uORs_file.close()
miss_counter = 0
correct_counter = 0

uORF_coords_file = open('D:/bioinformatic_institution/uORF_project/res_probs_2.txt', 'r')
uORF_coords_with_real_vars = open('D:/bioinformatic_institution/uORF_project/uORF_coords_with_real_vars_sins_clinsig_2.txt', 'w')

for each_str in uORF_coords_file:
    real_syn = 0
    real_mis = 0
    real_non = 0
    real_syn_sin = 0
    real_mis_sin = 0
    real_non_sin = 0
    clin_sig_syn = 0
    clin_sig_mis = 0
    clin_sig_non = 0
    inner_line = each_str.split(" ")
    coords = inner_line[0][3:] + "-" + inner_line[1] + "-" + inner_line[2]

    exons_intervals = exons_coords[coords]

    if coords in coords_dict and len(inner_line[5]) > 6:
        for each_var in coords_dict[coords]:
            if len(each_var[3]) == len(each_var[4]) == 1:
                for each_exon_interval in exons_intervals:
                    var_position = int(each_var[1]) - 1
                    if (each_exon_interval[0] <= var_position) and (var_position < each_exon_interval[1]):
                        rel_var_pos = each_exon_interval[2] + var_position - each_exon_interval[0]
                        if inner_line[4] == "-":
                            rel_var_pos = len(inner_line[5]) - 1 - rel_var_pos
                            each_var[3] = compl[each_var[3]]
                            each_var[4] = compl[each_var[4]]

                        try:
                            if inner_line[5][rel_var_pos] != each_var[3]:
                                miss_counter += 1
                                print(str(rel_var_pos) + ' ' + inner_line[5] + ' ' + inner_line[4] + ' ' + each_var[3] + ' ' + inner_line[5][rel_var_pos])
                            else:
                                correct_counter += 1
                        except:
                            continue
                        try:
                            if rel_var_pos % 3 == 0:
                                ref_codon = inner_line[5][rel_var_pos:rel_var_pos + 3]
                                ref_amino = gen_code[ref_codon]
                                var_codon = each_var[4] + inner_line[5][rel_var_pos + 1:rel_var_pos + 3]
                                var_amino = gen_code[var_codon]
                                if var_amino == "*":
                                    real_non += 1
                                    if int(each_var[16]) == 1:
                                        real_non_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_non += 1
                                elif ref_amino == var_amino:
                                    real_syn += 1
                                    if int(each_var[16]) == 1:
                                        real_syn_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_syn += 1
                                else:
                                    real_mis += 1
                                    if int(each_var[16]) == 1:
                                        real_mis_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_mis += 1
                            elif (rel_var_pos - 1) % 3 == 0:
                                ref_codon = inner_line[5][rel_var_pos - 1:rel_var_pos + 2]
                                ref_amino = gen_code[ref_codon]
                                var_codon = inner_line[5][rel_var_pos - 1] + each_var[4] + inner_line[5][rel_var_pos + 1]
                                var_amino = gen_code[var_codon]
                                if var_amino == "*":
                                    real_non += 1
                                    if int(each_var[16]) == 1:
                                        real_non_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_non += 1
                                elif ref_amino == var_amino:
                                    real_syn += 1
                                    if int(each_var[16]) == 1:
                                        real_syn_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_syn += 1
                                else:
                                    real_mis += 1
                                    if int(each_var[16]) == 1:
                                        real_mis_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_mis += 1
                            else:
                                ref_codon = inner_line[5][rel_var_pos - 2:rel_var_pos + 1]
                                ref_amino = gen_code[ref_codon]
                                var_codon = inner_line[5][rel_var_pos -2:rel_var_pos] + each_var[4]
                                var_amino = gen_code[var_codon]
                                if var_amino == "*":
                                    real_non += 1
                                    if int(each_var[16]) == 1:
                                        real_non_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_non += 1
                                elif ref_amino == var_amino:
                                    real_syn += 1
                                    if int(each_var[16]) == 1:
                                        real_syn_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_syn += 1
                                else:
                                    real_mis += 1
                                    if int(each_var[16]) == 1:
                                        real_mis_sin += 1
                                    if (each_var[13] == 'Pathogenic') or (each_var[13] == 'Likely pathogenic') or (each_var[13] == 'Pathogenic/Likely pathogenic'):
                                        clin_sig_mis += 1
                        except:
                            print(str(rel_var_pos) + ' ' + inner_line[5] + ' ' + inner_line[4])


        uORF_coords_with_real_vars.write(each_str[:-1] + " " + str(real_syn) + " " + str(real_mis) + " " + str(real_non) + " " + str(real_syn_sin) + " " + str(real_mis_sin) + " " + str(real_non_sin) + " " + str(clin_sig_syn) + " " + str(clin_sig_mis) + " " + str(clin_sig_non) + "\n")

print(miss_counter)
print(correct_counter)
uORF_coords_file.close()
uORF_coords_with_real_vars.close()



sequence = urllib.request.urlopen("http://togows.org/api/ucsc/hg38/chr16:57245264-57245294.fasta").read()
sequence = sequence.decode()
sequence = sequence.split('\n')
name = sequence[0]
whole_seq = "".join(sequence[1:])
print(name)
print(whole_seq)