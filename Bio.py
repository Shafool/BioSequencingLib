"""
Inicio clases
"""    

# Clase padre
class Secuencia:
# Constructor
    def __init__(self, seq, tipo_secuencia):
        self.seq = seq
        self.tipo_secuencia = tipo_secuencia

    def __len__(self):
        return len(self.seq)  
    
    def __str__(self):
        return self.seq
    
    def __getitem__(self,n):
        return  self.seq[n]

    def __getslice__(self, i, j):
        return self.seq[i:j]

# get y set
    def get_biotipo_seq(self):
        return self.tipo_secuencia

# Otras funciones

    # Cuenta el numero de veces que aparece una cadena en la secuencua
    def contar_ocurrencias(self, buscar_seq):
        return self.seq.count(buscar_seq)

    def imprimir_secuencia(self):
        print("Secuencia:", self.seq)

    def mostrar_info_secuencia(self):
        print("Secuencia:", self.seq, ", Biotipo:", self.tipo_secuencia)      

    # Recibe una secuencia y una lista de caracteres permitidos
    def validar_secuencia(self, secuencia, permitidos):
        seq = secuencia.upper()
        valido = True

        # Contamos todos los caracteres para ver si coincide con el largo de la secuencia
        for caracter in seq:
            if caracter not in permitidos:
                valido = False
                break
        return valido
    
    # Retorna un diccionario con la frecuencia de cada caracter en la secuencia
    def frecuencia(self):
        secuencia = self.seq
        dic = {}
        seq = secuencia.upper()

        for s in seq:
            if s in dic: dic[s] += 1
            else : dic[s] = 1
        
        return dic

    def invertir(self):
        secuencia = self.seq
        inv = ''
        for index in range(-1, -(len(secuencia) + 1), -1):
            inv += secuencia[index]
        return inv

    # Halla el complemento inverso de una secuencia, dado un diccionario, ejm ADN: {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
    def complemento_inverso(self, diccionario):
        secuencia = self.seq
        comp = ''
        for c in range(0, len(secuencia)):
            comp = diccionario[secuencia[c]] + comp
        return comp

# ========================================================

class Adn(Secuencia):
# Constructor personalizado
    # Validar la construccion del objeto
    def __init__(self, seq):
        t_seq = "ADN"
        secuencia = seq.upper()
        
        # Validar secuencia
        sec_valida = self.validar_adn(secuencia)
        if sec_valida == True:
            datos_validos = True
        else:
            datos_validos = False
            print("Error: uno o mas caracteres de la secuencia no coinciden con el tipo ADN")

        # Ingresar datos
        if datos_validos:
            self.tipo_secuencia = t_seq
            self.seq = secuencia

# get y set
    def set_biotipo_seq(self, bt):
        print("Error: el biotipo de esta secuencia debe ser solo de tipo ADN")

    def set_seq(self, seq):
        # valida que la nueva secuencia es adn antes de ingresarla
        secuencia = seq.upper()
        if validar_secuencia(secuencia, self.tipo_secuencia):
            self.seq = secuencia
        else:
            print("Error: esta secuencia debe ser de tipo ADN")

# Otras funciones
    def validar_adn(self, seq):
        nucleotidos = ['A', 'T', 'G', 'C']

        secuencia = seq.upper()
        valido = super().validar_secuencia(secuencia, nucleotidos)
        return valido
    
    # Retorna el porcentaje de Guanina y Citosina dada una secuencia de ADN, por defecto es self.seq
    def gc_porcentaje (self):
        adn_sec = self.seq
        gc_conteo = 0
        for s in adn_sec:
            if s in "GCgc": gc_conteo += 1
        return gc_conteo / len (adn_sec)

    # Retorna el porcentaje de CyG de subsecuencias no solapadas de tamaño k. Retorna una lista
    def gc_porcentaje_subsecuencia (self, k = 100):
        adn_sec = self.seq
        res = []

        for i in range (0, len(adn_sec) - k + 1, k):
            subseq = adn_sec[i:i+k]
            gc = gc_porcentaje(subseq)
            res.append(gc)
        return res

    # Complemento Inverso
    def complemento_inverso(self):
        adn_sec = self.seq
        assert self.validar_adn(adn_sec), "Error: Secuencia de ADN Inválida"

        conversion_dic = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
        comp = super().complemento_inverso(conversion_dic)
        return comp

    # Transcribe una secuencia de ADN a ARN
    def transcripcion(self):
        if (self.tipo_secuencia == "ADN"):
            return Arn(self.seq.replace("T","U"))
        else:
            return None

    # Traduce una secuencia de ADN directamente a una Proteina
    def traduccion(self, iniPos= 0):
        if (self.tipo_secuencia != "ADN"): return None
        seq_aa = ""
        for pos in range(iniPos,len(self.seq)-2,3):
            cod = self.seq[pos:pos+3]
            seq_aa += traducir_codon(cod)
        return Proteina(seq_aa)

    # Otras funciones extrañas xd
    # Frecuencia de cada codón codificando, dado como parámetro un aminoacido  
    def uso_codones(self, aa):
        dna_seq = (self.seq)
        seqm = dna_seq.upper()
        dic = {}
        total = 0
        for i in range(0, len(seqm)-2, 3):
            cod = seqm[i:i+3]
            if traducir_codon(cod) == aa:
                if cod in dic: 
                    dic[cod] += 1
                else: dic[cod] = 1
                total += 1
        if total >0:
            for k in dic:
                dic[k] /= total
        return dic


# Trabajo en proceso 

class Arn(Secuencia):
    # Constructor personalizado
    # Validar la construccion del objeto
    def __init__(self, seq):
        t_seq = "ARN"
        secuencia = seq.upper()
        
        # Validar secuencia
        sec_valida = self.validar_arn(secuencia)
        if sec_valida == True:
            datos_validos = True
        else:
            datos_validos = False
            print("Error: uno o mas caracteres de la secuencia no coinciden con el tipo ARN")

        # Ingresar datos
        if datos_validos:
            self.tipo_secuencia = t_seq
            self.seq = secuencia

# get y set
    def set_biotipo_seq(self, bt):
        print("Error: el biotipo de esta secuencia debe ser solo de tipo ARN")

    def set_seq(self, seq):
        # valida que la nueva secuencia es ARN antes de ingresarla
        secuencia = seq.upper()
        if validar_secuencia(secuencia, self.tipo_secuencia):
            self.seq = secuencia
        else:
            print("Error: esta secuencia debe ser de tipo ARN")

# Otras funciones
    def validar_arn(self, seq):
        nucleotidos = ['A', 'U', 'G', 'C']

        secuencia = seq.upper()
        valido = super().validar_secuencia(secuencia, nucleotidos)
        return valido
    
    # Hace la transcripción ARN -> ADN
    def transcripcion_inversa(self):
        arn_sec = self.seq
        assert self.validar_arn(arn_sec), "Error: Secuencia de ARN Inválida"
        return Adn(arn_sec.upper().replace("U","T"))

    # Traduce una secuencia de ARN a una Proteina
    def traduccion(self, iniPos= 0):
        if (self.tipo_secuencia != "ARN"): return None
        seq_aa = ""
        for pos in range(iniPos,len(self.seq)-2,3):
            cod = self.seq[pos:pos+3]
            # reemplazamos el uracilo con timina, pues la funcion traducir_codon solo funciona con adn
            cod = cod.replace("U","T")
            cod = cod.replace("u","T")
            cod = cod.upper()
            seq_aa += traducir_codon(cod)
        return Proteina(seq_aa)

class Proteina(Secuencia):
    # Constructor personalizado
    # Validar la construccion del objeto
    def __init__(self, seq):
        t_seq = "PROTEINA"
        secuencia = seq.upper()
        
        # Validar secuencia
        sec_valida = self.validar_prot(secuencia)
        if sec_valida == True:
            datos_validos = True
        else:
            datos_validos = False
            print("Error: uno o mas caracteres de la secuencia no coinciden con el tipo PROTEINA")

        # Ingresar datos
        if datos_validos:
            self.tipo_secuencia = t_seq
            self.seq = secuencia

# get y set
    def set_biotipo_seq(self, bt):
        print("Error: el biotipo de esta secuencia debe ser solo de tipo PROTEINA")

    def set_seq(self, seq):
        # valida que la nueva secuencia es PROTEINA antes de ingresarla
        secuencia = seq.upper()
        if validar_secuencia(secuencia, self.tipo_secuencia):
            self.seq = secuencia
        else:
            print("Error: esta secuencia debe ser de tipo PROTEINA")

# Otras funciones
    def validar_prot(self, seq):
        aminoacidos = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',]

        secuencia = seq.upper()
        valido = super().validar_secuencia(secuencia, aminoacidos)
        return valido
"""
Fin clases
"""       











"""
Sección de funciones utilitarias
"""        
# Funciones globales y utilitarias
def traducir_codon (cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", 
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None

def archivo2sec(filename):
    """ Reads a sequence from a multi-line text file. """
    fh = open(filename, "r")
    lines = fh.readlines()
    seq = ""
    for l in lines: 
        seq += l.replace("\n","")
    fh.close()
    return seq

def sec2archvo(seq, filename):
    """ Writes a sequence to file. """
    fh = open(filename, "w")
    fh.write(seq)
    fh.close()
    return None

"""Computes the six reading frames of a DNA sequence (including the reverse complement."""
def reading_frames (dna_seq):
    
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))    
    return res

def all_proteins_rf (aa_seq):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper() 
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_orfs (dna_seq):
    """Computes all possible proteins for all open reading frames."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames (dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots: res.append(p)
    return res

def all_orfs_ord (dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames (dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots: 
            if len(p) > minsize: insert_prot_ord(p, res)
    return res

# Daniel

def find_pattern_re (seq, pat):
    from re import search
    mo = search(pat, seq)
    if (mo != None):
        return mo.span()[0]
    else: 
        return -1

def find_all_occurrences_re (seq, pat):
    from re import finditer
    mos = finditer(pat, seq)
    res = []
    for x in mos:
        res.append(x.span()[0])
    return res

def find_all_overlap(seq, pat):
    return find_all_occurrences_re(seq, "(?="+pat+")")

def validate_dna_re (seq):
    from re import search
    if search("[^ACTGactg]", seq) != None:
        return False
    else :
        return True

def translate_codon_re (cod):
    import re
    if re.search("GC.", cod): aa = "A"
    elif re.search("TG[TC]", cod): aa = "C"
    elif re.search("GA[TC]", cod): aa = "D"
    elif re.search("GA[AG]", cod): aa = "E"
    elif re.search("TT[TC]", cod): aa = "F"
    elif re.search("GG.", cod): aa = "G"
    elif re.search("CA[TC]", cod): aa = "H"
    elif re.search("AT[TCA]", cod): aa = "I"
    elif re.search("AA[AG]", cod): aa = "K"
    elif re.search("TT[AG]|CT.", cod): aa = "L"
    elif re.search("ATG", cod): aa = "M"
    elif re.search("AA[TC]", cod): aa = "N"
    elif re.search("CC.", cod): aa = "P"
    elif re.search("CA[AG]", cod): aa = "Q"
    elif re.search("CG.|AG[AG]", cod): aa = "R"
    elif re.search("TC.|AG[TC]", cod): aa = "S"
    elif re.search("AC.", cod): aa = "T"
    elif re.search("GT.", cod): aa = "V"
    elif re.search("TGG", cod): aa = "W"
    elif re.search("TA[TC]", cod): aa = "Y"
    elif re.search("TA[AG]|TGA", cod): aa = "_";
    else : aa = ""
    return aa

def largest_protein_re (seq_prot):
    import re
    mos = re.finditer("M[^_]∗_", seq_prot)
    sizem = 0
    lprot = ""
    for x in mos:
        ini = x.span()[0]
        fin = x.span()[1]
        s = fin - ini + 1
        if s > sizem:
            lprot = x.group()
            sizem = s
        return lprot

def find_zync_finger(seq):
    from re import search
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else :
        return -1

def find_prosite(seq, profile):
    from re import search
    regexp = profile.replace("−","")
    regexp = regexp.replace("x",".")
    regexp = regexp.replace("(","{")
    regexp = regexp.replace(")","}")
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else :
        return -1
def cut_positions (enzyme, sequence):
    from re import finditer
    cutpos = enzyme.find("^")
    regexp = iub_to_RE(enzyme)
    matches = finditer(regexp, sequence)
    locs = [ ]
    for m in matches:
        locs.append(m.start() + cutpos)
    return locs

def cut_subsequences (locs, sequence):
    res = []
    positions = locs
    positions.insert(0,0)
    positions.append( len (sequence))
    for i in range( len (positions)-1):
        res.append(sequence[positions[i]:positions[i+1]])
    return res

def iub_to_RE (iub):
    dic = {"A":"A", "C":"C", "G":"G", "T":"T", "R":"[GA]", "Y":"[CT]", "M":"[AC]", "K":"[GT]", "S":"[GC]", "W": "[AT]", "B":"[CGT]", "D":"[AGT]", "H":"[ACT]", "V":"[ACG]", "N":"[ACGT]"}
    site = iub.replace("^","")
    regexp = ""
    for c in site:
        regexp += dic[c]
    return regexp

def search_first_occ(seq, pattern):
    found = False
    i=0

    while i <= len (seq)-len (pattern) and not found:
        j=0
        while j < len (pattern) and pattern[j]==seq[i+j]:
            j=j+1
        if j== len (pattern): found = True
        else : i += 1
    if found: return i
    else:return -1
    
    
def search_all_occurrences(seq, pattern):
    res = []
    for i in range( len (seq)-len (pattern) +1):
        j=0
        while j < len (pattern) and pattern[j]==seq[i+j]:
            j=j+1
        if j == len (pattern):
            res.append(i)
    return res


class BoyerMoore:
    def __init__( self , alphabet, pattern):
        self .alphabet = alphabet
        self .pattern = pattern
        self .preprocess()
    def preprocess( self ):
        self .process_bcr()
        self .process_gsr()
    def process_bcr( self ):
        self .occ = {}
        for symb in self .alphabet:
            self.occ[symb] = -1
        for j in range( len ( self .pattern)):
            c = self .pattern[j]
            self .occ[c] = j
    def process_gsr( self ):
        self .f = [0] * ( len ( self .pattern)+1)
        self .s = [0] * ( len ( self .pattern)+1)
        i = len ( self .pattern)
        j = len ( self .pattern)+1
        self .f[i] = j
        while i>0:
            while j<= len ( self .pattern) and self .pattern[i-1] != self.pattern[j-1]:
                if self .s[j] == 0: self .s[j] = j-i;
                j = self .f[j]
            i -= 1
            j -= 1
            self .f[i] = j
        j = self .f[0]
        for i in range( len ( self .pattern)):
            if self .s[i] == 0: self .s[i] = j
            if i == j: j = self .f[j]
    def search_pattern( self , text):
        res = []
        i=0
        while i <= len (text) - len ( self .pattern):
            j= len ( self .pattern)- 1
            while j>=0 and self .pattern[j]==text[j+i]: j -= 1
            if (j<0):
                res.append(i)
                i += self .s[0]
            else :
                c = text[j+i]
                i += max( self .s[j+1], j- self .occ[c])
        return res

def overlap(s1, s2):
    maxov = min( len (s1), len (s2))
    for i in range(maxov,0,-1):
        if s1[-i:] == s2[:i]: return i
    return 0
class Automata:
    def __init__( self , alphabet, pattern):
        self .numstates = len (pattern) + 1
        self .alphabet = alphabet
        self .transition_table = {}
        self .build_transition_table(pattern)
    def build_transition_table( self , pattern):
        for q in range( self .numstates):
            for a in self .alphabet:
                prefix = pattern[0:q] + a
                self .transition_table[(q,a)] = overlap (prefix, pattern)
    def print_automata( self ):
        print ("States: " , self .numstates)
        print ("Alphabet: " , self .alphabet)
        print ("Transition table:")
        for k in self .transition_table.keys():
            print (k[0], ",", k[1], " −> ", self .transition_table[k])
            
    def next_state( self , current, symbol):
        return self .transition_table.get((current, symbol))
    
    def apply_seq( self , seq):
        q=0
        res = [q]
        for c in seq:
            q = self .next_state(q, c)
            res.append(q)
        return res

    def occurrences_pattern( self , text):
        q=0
        res = []
        for i in range( len (text)):
            q = self .next_state(q, text[i])
            if q == self .numstates-1:
                res.append(i - self .numstates + 2)
        return res

# ==================================================================

