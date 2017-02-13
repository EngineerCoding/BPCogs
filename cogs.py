import os

def OpenFile(of_file):
    with open(of_file,'r') as f:
        inhoud = f.readlines()
        f.close()
        index = 0
        while index < len(inhoud):
            inhoud[index] = inhoud[index].rstrip("\n")
            index += 1
        return inhoud

def Database():
    database_code = OpenFile("Database_commando.txt")
    
    for x in database_code:
        print("Database", x, "in aanmaak.")
        os.system(x)

def Blast(organisme_1, organisme_2):
    blastnamen = []
    count = 1
    for x in organisme_1:
        for z in organisme_2:
            if x == z:
                pass
            else:
                #print("Blast", count, "van de 90: \t", x, "tegen", z)
                count += 1
                #os.system("blastp -outfmt 6 -query %s -db %s -evalue 0.0000000001 -num_threads 3 | awk '{print $1,$2}' | sed 's/ /;/g' > %s.txt"%(x, z, x+"_tegen_"+z))
                blastnamen.append(x+"_tegen_"+z+".txt")
    return(blastnamen)


def Original(blast_resultaten):
    bestanden = []
    count = 1
    for x in blast_resultaten:
        index = 0
        nieuwe_naam_deel1 = x.split("_tegen_")[0]
        nieuwe_naam_deel2 = x.split("_tegen_")[1]
        nieuwe_naam = nieuwe_naam_deel1 + "_best_tegen_" + nieuwe_naam_deel2
        bestanden.append(nieuwe_naam)
##        hit = OpenFile(x)
##        best_hit = []
##        vorige_hit = ["deze index positie wordt vervangen"]
##        for x in hit:
##            naam = x.split(";")[0]
##            if naam != vorige_hit[0]:
##                best_hit.append(hit[index])
##                vorige_hit.pop(0)
##                vorige_hit.append(naam)
##            index += 1
##        Make_file(nieuwe_naam, best_hit)
    return(bestanden)


def Make_file(name, inhoud):
    outfile = open(name,'a')
    for x in inhoud:
        outfile.write(x + "\n")
    outfile.close

def Reverse(omdraaien):
    omgekeerd = []
    for z in omdraaien:
        hit1 = z.split(";")[0]
        hit2 = z.split(";")[1]
        omgekeerd.append(hit2 + ";" + hit1)
    return(omgekeerd)


def Bi(bestanden):
    verbruikte_bestanden = []
    for bestand in bestanden:
        split1 = bestand.split("_best_tegen_")[1]
        split2 = split1.split(".txt")[0]
        split3 = bestand.split("_best_tegen_")[0]
        bestand2 = split2 + "_best_tegen_" + split3 + ".txt"
        if bestand in verbruikte_bestanden or bestand2 in verbruikte_bestanden:
            pass
        else:
            print(bestand)
            verbruikte_bestanden.append(bestand)
            verbruikte_bestanden.append(bestand2)
            A = OpenFile(bestand)
            B = OpenFile(bestand2)
            B_omgekeerd = Reverse(B)
            bidirectional = []
            bidirectional.append("> bidirectionalhits tussen " + bestand)
            for x in A:
                for z in B_omgekeerd:
                    if x == z:
                        bidirectional.append(x)
            Make_file("bidirectionalhits.txt", bidirectional) 
            
  

def Main():
    organismen = OpenFile("Organismen.txt")
    #Database()
    blast_resultaten = Blast(organismen, organismen)
    bestanden = Original(blast_resultaten)
    Bi(bestanden)
    

Main()
