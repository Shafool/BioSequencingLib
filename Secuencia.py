#Clase padre
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
    
    # Hace la transcripción ADN -> ARN y retorna un string
    def transcripcion (self):
        adn_sec = self.seq
        assert self.validar_adn(adn_sec), "Error: Secuencia de ADN Inválida"
        return adn_sec.upper().replace("T","U")

    # Complemento Inverso
    def complemento_inverso(self):
        adn_sec = self.seq
        assert self.validar_adn(adn_sec), "Error: Secuencia de ADN Inválida"

        conversion_dic = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
        comp = super().complemento_inverso(conversion_dic)
        return comp

# Trabajo en proceso
class Arn(Secuencia):

class Proteina(Secuencia):