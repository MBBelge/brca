import matplotlib.pyplot as plt
import re

def varyant_bul(substr, infile, outfile):
    ekzonik_count = intronik_count = deletion_count = insertion_count = duplication_count = others_count = at_count = ta_count = gt_count = tg_count = gc_count = cg_count = ac_count = ca_count = ag_count = ga_count = ct_count = tc_count = 0
    set = 0
    exonsozluk = {}
    onceki = ''

    with open(infile) as a, open(outfile, 'w') as b:
        for line in a:
            elements = line.strip().split('\t')
            if elements[1] == substr and elements[3].startswith('Pathogenic') and elements[3] != 'Likely pathogenic':
                variant_name = elements[0].split(' ')
                var_name = variant_name[0].split('.')
                var_name = var_name[2]
                var_name = re.sub(r'\d*', '',var_name)
                var_name = re.sub(r'\**', '', var_name)
                var_name = re.sub(r'\_*', '', var_name)
                var_name = re.sub(r'\-*', '', var_name)
                if ";" not in var_name:
                    if "del" in var_name:
                        deletion_count += 1
                    if "ins" in var_name:
                        insertion_count += 1
                    if "dup" in var_name:
                        duplication_count += 1
                    if ">" in var_name:
                        others_count += 1
                        if "A>T" in var_name:
                            at_count += 1
                        if "T>A" in var_name:
                            ta_count += 1
                        if "G>T" in var_name:
                            gt_count += 1
                        if "T>G" in var_name:
                            tg_count += 1
                        if "G>C" in var_name:
                            gc_count += 1
                        if "C>G" in var_name:
                            cg_count += 1
                        if "A>C" in var_name:
                            ac_count += 1
                        if "C>A" in var_name:
                            ca_count += 1
                        if "A>G" in var_name:
                            ag_count += 1
                        if "G>A" in var_name:
                            ga_count += 1
                        if "C>T" in var_name:
                            ct_count += 1
                        if "T>C" in var_name:
                            tc_count += 1
                    if len(variant_name) == 2:
                        ekzonik_count = ekzonik_count + 1
                        ekzonik_sozluk = {elements[0]: [elements[5], elements[6]]}
                        b.write(('\t'.join([elements[0], elements[5], elements[6]]) + '\n'))
                        print(ekzonik_sozluk)
                    elif len(variant_name) == 1:
                        intronik_count = intronik_count + 1
                        sozluk = {}
                        sozluk = {elements[0]: [elements[5], elements[6]]}

        piechart(deletion_count, insertion_count, duplication_count,others_count)
        piechart2(at_count, ta_count, gt_count, tg_count, gc_count, cg_count, ac_count, ca_count, ag_count, ga_count,
             ct_count, tc_count)
        piechart3(intronik_count,ekzonik_count)

    print(" ")
    with open('C:\\Users\\MBB\\Desktop\\proj\\exon.txt', 'r') as c:
        with open('C:\\Users\\MBB\\Desktop\\proj\\output.txt', 'r') as d:
            with open('C:\\Users\\MBB\\Desktop\\proj\\exonVariant.txt', 'w') as e:
                with open('C:\\Users\\MBB\\Desktop\\proj\\exonSayac.txt', 'w') as f:
                    for line in c:
                        elements = line.strip().split('\t')
                        xkoor = elements[9].strip().split(",")
                        ykoor = elements[10].strip().split(",")
                        for line2 in d:
                            element = line2.strip().split('\t')
                            koor = re.match(r'\d+', element[2])
                            loop = 0
                            for i in xkoor:
                                i = int(i)
                                j = int(ykoor[loop])
                                if onceki == '':
                                    onceki = str(i)
                                if (int(koor.group(0)))>i and (int(koor.group(0)))<j:
                                    e.write('\t'.join([element[0], element[1], element[2]]) + '\n')
                                    print('Success')
                                    if i in exonsozluk.keys():
                                        count = exonsozluk.get(i)
                                        exonsozluk[i] = count+1
                                    else:
                                        exonsozluk = {i:1}
                                    break
                                else:
                                    print('Failure')
                                if str(exonsozluk.get(i)) != 'None':
                                    f.write(str(i) + ' ve ' + str(j) + ' Arasinda ' ' : ' + str(exonsozluk.get(i)) + '\n')
                                loop += 1
            ekzon_list = line.strip().split('\t')
            plt.scatter(xkoor, ykoor, alpha=1, s=1)
            plt.show()
            print(ekzon_list)

    with open('C:\\Users\\MBB\\Desktop\\proj\\count.txt', 'w') as c:
            c.write("Eksonik Sayisi : {0}\nIntronik Sayisi : {1}\nDeletion Sayisi : {2}\nInsertion Sayisi : {3}\nDuplication Sayisi : {4}\nSubstitute Sayisi : {5}\nA>T Sayisi : {6}\nT>A Sayisi : {7}\nG>T Sayisi : {8}\nT>G Sayisi : {9}\nG>C Sayisi : {10}\nC>G Sayisi : {11}\nA>C Sayisi : {12}\nC>A Sayisi : {13}\nA>G Sayisi : {14}\nG>A Sayisi : {15}\nC>T Sayisi : {16}\nT>C Sayisi : {17}".format(ekzonik_count,intronik_count,deletion_count,insertion_count,duplication_count,others_count,at_count, ta_count, gt_count, tg_count, gc_count, cg_count, ac_count, ca_count, ag_count, ga_count,ct_count, tc_count))

def piechart(deletion_count,insertion_count,duplication_count,others_count):
    labels = 'Deletion', 'Insertion', 'Duplication', 'Substitions'
    sizes = [deletion_count, insertion_count, duplication_count,others_count]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
    ax1.axis('equal')
    plt.show()

def piechart3(intronik_count,ekzonik_count):
    labels = 'Intronik Sayisi','Ekzonik Sayisi'
    sizes = [intronik_count,ekzonik_count]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
    ax1.axis('equal')
    plt.show()

def piechart2(at_count, ta_count, gt_count, tg_count, gc_count, cg_count, ac_count, ca_count, ag_count, ga_count, ct_count, tc_count):
    labels = 'A>T', 'T>A', 'G>T', 'T>G', 'G>C', 'C>G', 'A>C', 'C>A', 'A>G', 'G>A', 'C>T', 'T>C'
    sizes = [at_count, ta_count, gt_count, tg_count, gc_count, cg_count, ac_count, ca_count, ag_count, ga_count, ct_count, tc_count]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
    ax1.axis('equal')
    plt.show()

varyant_bul('BRCA1', 'C:\\Users\MBB\\Desktop\\proj\\clinvar_result.txt', 'C:\\Users\\MBB\\Desktop\\proj\\output.txt')
