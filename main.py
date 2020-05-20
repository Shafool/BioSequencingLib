import Bio
from Bio import Secuencia
from Bio import Adn
from Bio import Arn

# Especificando la secuencia
sec_input = input("> Ingrese una secuencia de ADN: ")
adnseq = Adn(sec_input)  

# Llamando y probando metodos ADN
adnseq.mostrar_info_secuencia()
print("\tComplemento inverso:", adnseq.complemento_inverso())
print("\tPorcentaje de C y G:",adnseq.gc_porcentaje())
print("\tTraduccion (Proteina):",adnseq.traduccion().seq)
print("\tTranscripcion (ADN -> ARN):",adnseq.transcripcion().seq)

print("")

# Especificando la secuencia ARN
sec_input = input("> Ingrese una secuencia de ARN: ")
arnseq = Arn(sec_input)  

# Llamando y probando metodos ADN
arnseq.mostrar_info_secuencia()
print("\tTraduccion (Proteina):",arnseq.traduccion().seq)
print("\tTranscripcion inversa (ARN -> ADN):",arnseq.transcripcion_inversa().seq)