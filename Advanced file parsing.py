#Thao My Nguyen
#Advanced file parsing
#Email: mnguye39@uncc.edu


aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 
'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'], 
'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 
'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'], 
'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}


def nuc_seq_dict(file):
    fh = open("Mdomestica_491_v1.1.cds_primaryTranscriptOnly.txt","r")

    fn=fh.readlines()


    seq_dct={}                                                  
    count=1                                                     
    for line in fn:

        if line.startswith('>MD'):
            header=line.rstrip('\n').strip(">MD")
            sequence=(fn[count].strip("\n"))

            seq_dct[header]=sequence
            count+=2                                          
        else:
            continue
   
    return seq_dct

#function to count total number of amino acid 
def aa_count(a):
    #set amino acid to 0 then doing for loop to count number of each amino acid 
    Leu=0
    Met=0
    Phe=0
    Cys=0
    Tyr=0
    Trp=0
    Pro=0
    His=0
    Gln=0
    Arg=0
    Ile=0
    Thr=0
    Asn=0
    Lys=0
    Ser=0
    Val=0
    Ala=0
    Asp=0
    Glu=0
    Gly=0
    Asterick=0
    seq=a

    for key,value in seq.items():


	    for i in range (0,len(value),3):
	        print(value[i:i+3])

	        if value[i:i+3]=="CTT" or'TTA' or 'TTG'or 'CTC' or 'CTA' or 'CTG':
	            Leu+=1
	        print("Amino acid Leu count:", Leu)
	        if value[i:i+3]=="ATG":
	            Met+=1
	        print("Amino acid Met count:", Met)
	        if value[i:i+3]=='TTT'or 'TTC':
	            Phe+=1
	        print("Amino acid Phe count:", Phe)
	        if value[i:i+3]=='TGT' or 'TGC':
	            Cys+=1
	        print("Amino acid Cys count:", Cys)
	        if value[i:i+3]=='TAC' or 'TAT':
	            Tyr+=1
	        print("Amino acid Tyr count:", Tyr)
	        if value[i:i+3]=='TGG':
	            Trp+=1
	        print("Amino acid Trp count:", Trp)
	        if value[i:i+3]=='CCT' or 'CCC' or 'CCA' or 'CCG':
	            Pro+=1
	        print("Amino acid Pro count:", Pro)
	        if value[i:i+3]=='CAT' or 'CAC':
	            His+=1
	        print("Amino acid His count:", His)
	        if value[i:i+3]=='CAA' or'CAG':
	            Gln+=1
	        print("Amino acid Gln count:", Gln)
	        if value[i:i+3]=='CGT' or 'CGC' or 'CGA' or 'CGG' or 'AGA' or 'AGG':
	            Arg+=1
	        print("Amino acid Arg count:", Arg)
	        if value[i:i+3]=='ATT' or'ATC' or 'ATA':
	            Ile+=1
	        print("Amino acid Ile count:", Ile)
	        if value[i:i+3]=='ACT' or 'ACC' or 'ACA' or 'ACG':
	            Thr+=1
	        print("Amino acid Thr count:", Thr)
	        if value[i:i+3]=='AAT' or'AAC':
	            Asn+=1
	        print("Amino acid Asn count:", Asn)
	        if value[i:i+3]=='AAA' or'AAG':
	            Lys+=1
	        print("Amino acid Lys count:", Lys)
	        if value[i:i+3]=='AGT' or 'AGC' or'TCT' or 'TCC' or 'TCA' or 'TCG':
	            Ser+=1
	        print("Amino acid Ser count:", Ser)
	        if value[i:i+3]=='GTT' or 'GTC' or 'GTA' or'GTG':
	            Val+=1
	        print("Amino acid Val count:", Val)
	        if value[i:i+3]=='GCT' or 'GCC' or 'GCA' or 'GCG':
	            Ala+=1
	        print("Amino acid Ala count:", Ala)
	        if value[i:i+3]=='GAT' or 'GAC':
	            Asp+=1
	        print("Amino acid Asp count:", Asp)
	        if value[i:i+3]=='GAA'or 'GAG':
	            Glu+=1
	        print("Amino acid Glu count:", Glu)
	        if value[i:i+3]=='GGT' or 'GGC' or 'GGA' or 'GGG':
	            Gly+=1
	        print("Amino acid Gly count:", Gly)
	        if value[i:i+3]=='TAA' or 'TAG' or'TGA':
	            Asterick+=1
	        print("Amino acid Asterick count:", Asterick)

def codon_count(a):
	#set dictionary for total codon count
	codon_count= {'ATG':0, 'TTT':0, 'TTC':0,'TTA':0, 'TTG':0, 'CTT':0, 'CTC':0, 'CTA':0, 'CTG':0,'TGT':0, 'TGC':0,'TAC':0, 'TAT':0,'TGG':0, 'CCT':0, 'CCC':0, 'CCA':0, 'CCG':0, 'CAT':0, 'CAC':0, 
'CAA':0, 'CAG':0, 'CGT':0, 'CGC':0, 'CGA':0, 'CGG':0, 'AGA':0, 'AGG':0, 'ATT':0, 'ATC':0, 'ATA':0, 'ACT':0, 'ACC':0, 'ACA':0, 'ACG':0, 
'AAT':0, 'AAC':0, 'AAA':0, 'AAG':0, 'AGT':0, 'AGC':0, 'TCT':0, 'TCC':0, 'TCA':0, 'TCG':0, 'GTT':0, 'GTC':0, 'GTA':0, 'GTG':0, 
'GCT':0, 'GCC':0, 'GCA':0, 'GCG':0, 'GAT':0, 'GAC':0, 'GAA':0, 'GAG':0, 'GGT':0, 'GGC':0, 'GGA':0, 'GGG':0, 'TAA':0,'TAG':0,'TGA':0}
	print(codon_count)

	for key,value in seq.items():


		for i in range (0,len(value),3):
			
			print(value[i:i+3])
			if value[i:i+3]=="ATG":
				codon_count['ATG']+=1
			print("ATG:", codon_count['ATG'])
			if value[i:i+3]=='TTT':
				codon_count['TTT']+=1
			print("TTT:", codon_count['TTT'])
			if value[i:i+3]=="TTC":
				codon_count['TTC']+=1
			print("TTC:", codon_count['TTC'])
	        
			if value[i:i+3]=='TTA':
				codon_count['TTA']+=1
			print("TTA:", codon_count['TTA'])
			if value[i:i+3]=="TTG":
				codon_count['TTG']+=1
			print("TTG:", codon_count['TTG'])
			if value[i:i+3]=='CTT':
				codon_count['CTT']+=1
			print("CTT:", codon_count['CTT'])
			if value[i:i+3]=="CTC":
				codon_count['CTC']+=1
			print("CTC:", codon_count['CTC'])
			if value[i:i+3]=='CTA':
				codon_count['CTA']+=1
			print("CTA:", codon_count['CTA'])
			if value[i:i+3]=="CTG":
				codon_count['CTG']+=1
			print("CTG:", codon_count['CTG'])
	        
			if value[i:i+3]=='TGT':
				codon_count['TGT']+=1
			print("TGT:", codon_count['TGT'])
			if value[i:i+3]=="TGC":
				codon_count['TGC']+=1
			print("TGC:", codon_count['TGC'])
	        
			if valvalueue[i:i+3]=='TAC':
				codon_count['TAC']+=1
			print("TAC:", codon_count['TAC'])
			if value[i:i+3]=="TAT":
				codon_count['TAT']+=1
			print("TAT:", codon_count['TAT'])
	        
			if value[i:i+3]=='TGG':
				codon_count['TGG']+=1
			print("TGG:", codon_count['TGG'])
	        
			if value[i:i+3]=="CCT":
				codon_count['CCT']+=1
			print("CCT:", codon_count['CCT'])
			if value[i:i+3]=='CCC':
				codon_count['CCC']+=1
			print("CCC:", codon_count['CCC'])
			if value[i:i+3]=="CCA":
				codon_count['CCA']+=1
			print("CCA:", codon_count['CCA'])
			if value[i:i+3]=='CCG':
				codon_count['CCG']+=1
			print("CCG:", codon_count['CCG'])
	        
			if value[i:i+3]=="CAT":
				codon_count['CAT']+=1
			print("CAT:", codon_count['CAT'])
			if value[i:i+3]=='CAC':
				codon_count['CAC']+=1
			print("CAC:", codon_count['CAC'])
	        
			if value[i:i+3]=="CAA":
				codon_count['CAA']+=1
			print("CAA:", codon_count['CAA'])
			if value[i:i+3]=='CAG':
				codon_count['CAG']+=1
			print("CAG:", codon_count['CAG'])
	        
			if value[i:i+3]=="CGT":
				codon_count['CGT']+=1
			print("CGT:", codon_count['CGT'])
			if value[i:i+3]=='CGC':
				codon_count['CGC']+=1
			print("CGC:", codon_count['CGC'])
			if value[i:i+3]=="CGA":
				codon_count['CGA']+=1
			print("CGA:", codon_count['CGA'])
			if value[i:i+3]=='CGG':
				codon_count['CGG']+=1
			print("CGG:", codon_count['CGG'])
			if value[i:i+3]=="AGA":
				codon_count['AGA']+=1
			print("AGA:", codon_count['AGA'])
			if value[i:i+3]=='AGG':
				codon_count['AGG']+=1
			print("AGG:", codon_count['AGG'])
	        
			if value[i:i+3]=="ATT":
				codon_count['ATT']+=1
			print("ATT:", codon_count['ATT'])
			if value[i:i+3]=='ATC':
				codon_count['ATC']+=1
			print("ATC:", codon_count['ATC'])
			if value[i:i+3]=="ATA":
				codon_count['ATA']+=1
			print("ATA:", codon_count['ATA'])
	        
			if value[i:i+3]=='ACT':
				codon_count['ACT']+=1
			print("ACT:", codon_count['ACT'])
			if value[i:i+3]=="ACC":
				codon_count['ACC']+=1
			print("ACC:", codon_count['ACC'])
			if value[i:i+3]=='ACA':
				codon_count['ACA']+=1
			print("ACA:", codon_count['ACA'])
			if value[i:i+3]=="ACG":
				codon_count['ACG']+=1
			print("ACG:", codon_count['ACG'])
	        
			if value[i:i+3]=='AAT':
				codon_count['AAT']+=1
			print("AAT:", codon_count['AAT'])
			if value[i:i+3]=="AAC":
				codon_count['AAC']+=1
			print("AAC:", codon_count['AAC'])
	        
			if value[i:i+3]=='AAA':
				codon_count['AAA']+=1
			print("AAA:", codon_count['AAA'])
			if value[i:i+3]=="AAG":
				codon_count['AAG']+=1
			print("AAG:", codon_count['AAG'])
	        
			if value[i:i+3]=='AGT':
				codon_count['AGT']+=1
			print("AGT:", codon_count['AGT'])
			if value[i:i+3]=="AGC":
				codon_count['AGC']+=1
			print("AGC:", codon_count['AGC'])
			if value[i:i+3]=='TCT':
				codon_count['TCT']+=1
			print("TCT:", codon_count['TCT'])
			if value[i:i+3]=="TCC":
				codon_count['TCC']+=1
			print("TCC:", codon_count['TCC'])
			if value[i:i+3]=='TCA':
				codon_count['TCA']+=1
			print("TCA:", codon_count['TCA'])
			if value[i:i+3]=="TCG":
				codon_count['TCG']+=1
			print("TCG:", codon_count['TCG'])

			if value[i:i+3]=='GTT':
				codon_count['GTT']+=1
			print("GTT:", codon_count['GTT'])
			if value[i:i+3]=="GTC":
				codon_count['GTC']+=1
			print("GTC:", codon_count['GTC'])
			if value[i:i+3]=='GTA':
				codon_count['GTA']+=1
			print("GTA:", codon_count['GTA'])
			if value[i:i+3]=='GTG':
				codon_count['GTG']+=1
			print("GTG:", codon_count['GTG'])

			if value[i:i+3]=="GCT":
				codon_count['GCT']+=1
			print("GCT:", codon_count['GCT'])
			if value[i:i+3]=='GCC':
				codon_count['GCC']+=1
			print("GCC:", codon_count['GCC'])
			if value[i:i+3]=="GCA":
				codon_count['GCA']+=1
			print("GCA:", codon_count['GCA'])
			if value[i:i+3]=='GCG':
				codon_count['GCG']+=1
			print("GCG:", codon_count['GCG'])

			if value[i:i+3]=="GAT":
				codon_count['GAT']+=1
			print("GAT:", codon_count['GAT'])
			if value[i:i+3]=='GAC':
				codon_count['GAC']+=1
			print("GAC:", codon_count['GAC'])

			if value[i:i+3]=="GAA":
				codon_count['GAA']+=1
			print("GAA:", codon_count['GAA'])
			if value[i:i+3]=='GAG':
				codon_count['GAG']+=1
			print("GAG:", codon_count['GAG'])

			if value[i:i+3]=="GGT":
				codon_count['GGT']+=1
			print("GGT:", codon_count['GGT'])
			if value[i:i+3]=='GGC':
				codon_count['GGC']+=1
			print("GGC:", codon_count['GGC'])
			if value[i:i+3]=="GGA":
				codon_count['GGA']+=1
			print("GGA:", codon_count['GGA'])
			if value[i:i+3]=='GGG':
				codon_count['GGG']+=1
			print("GGG:", codon_count['GGG'])

			if value[i:i+3]=="TAA":
				codon_count['TAA']+=1
			print("TAA:", codon_count['TAA'])
			if value[i:i+3]=='TAG':
				codon_count['TAG']+=1
			print("TAG:", codon_count['TAG'])
			if value[i:i+3]=="TGA":
				codon_count['TGA']+=1
			print("TGA:", codon_count['TGA'])

		print(codon_count.items())

def main():

	a=nuc_seq_dict("Mdomestica_491_v1.1.cds_primaryTranscriptOnl.txt")
	aa_count(a)
	codon_count(a)


if __name__ == "__main__":
	main()