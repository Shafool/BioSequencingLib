from Secuencia import Secuencia
from Secuencia import Adn

# Especificando la secuencia
sec_input = input("Ingrese una secuencia: ")
s1 = Adn(sec_input)  

# Llamando y probando metodos
s1.imprimir_secuencia()
s1.mostrar_info_secuencia()
print("Complemento inverso:", s1.complemento_inverso())
print("Transcripcion:", s1.transcripcion())
print("Porcentaje de C y G:",s1.gc_porcentaje())