# Formatos SAM, BAM y CRAM

El dia de hoy obtendremos información esencial sobre los formatos SAM, BAM y CRAM utilizados en bioinformática para almacenar y procesar datos de alineamientos de secuencias. 

---

### ¿Qué es un archivo SAM?

El formato **SAM** (Sequence Alignment/Map) es un formato de texto delimitado por tabulaciones que contiene información sobre los alineamientos. Consiste en:

1. Una sección de encabezado, que incluye metadatos.
2. Una sección de alineamientos, donde cada línea proporciona información detallada de un alineamiento.

**Referencias útiles:**
- [Especificación oficial del formato SAM](http://samtools.github.io/hts-specs/SAMv1.pdf)
- [Especificación de etiquetas SAM](http://samtools.github.io/hts-specs/SAMtags.pdf)

---

### ¿Qué es un archivo BAM?

El formato **BAM** es la versión binaria y comprimida del archivo SAM. Generalmente, los archivos BAM están ordenados por coordenadas o nombres de lectura.

- **Ordenados por coordenadas:** Permite consultas rápidas basadas en ubicaciones.
- **Ordenados por nombres:** Facilita el acceso a pares de lectura con el mismo nombre.

---

### ¿Qué es un archivo CRAM?

El formato **CRAM** es una versión más eficiente de BAM. Requiere almacenar la referencia genómica por separado, lo que reduce significativamente el tamaño del archivo. Sin embargo, esto añade complejidad ya que la referencia debe estar disponible para decodificar el archivo.

**Ejemplo de creación:**
```bash
bwa mem $REF $R1 $R2 | samtools sort --reference $REF -O CRAM > bwa.cram
```

---

### Uso de los archivos SAM

Los archivos SAM están diseñados para:

1. **Almacenar alineamientos de manera estándar y eficiente.**
2. **Permitir acceso rápido por coordenadas.**

Por ejemplo, se puede recuperar alineamientos que se superpongan con una coordenada específica en milisegundos usando un archivo BAM ordenado.

---

### Creación de archivos SAM/BAM/CRAM

**Pasos básicos:**

1. Crear un archivo SAM:
   ```bash
   bwa mem $REF $R1 $R2 > alignments.sam
   ```
2. Convertir SAM a BAM:
   ```bash
   samtools sort alignments.sam > alignments.bam
   ```
3. Indexar el archivo BAM:
   ```bash
   samtools index alignments.bam
   ```

**Crear un archivo CRAM:**
```bash
bwa mem $REF $R1 $R2 | samtools sort --reference $REF -O cram > bwa.cram
samtools index bwa.cram
```
---

### Ejemplo práctico: Creación de un archivo BAM

A continuación, se muestra un ejemplo práctico para crear un archivo BAM usando datos públicos:

1. **Definir los parámetros principales:**
   ```bash
   ACC=AF086833  # Número de acceso del genoma de referencia
   SRR=SRR1972739  # Número de acceso del conjunto de datos
   N=10000  # Número de lecturas a extraer
   ```

2. **Preparar el genoma de referencia:**
   ```bash
   mkdir -p refs
   efetch -db nucleotide -format gb -id $ACC > refs/$ACC.gb
   cat refs/$ACC.gb | seqret -filter -feature -osformat fasta > refs/$ACC.fa
   samtools faidx refs/$ACC.fa
   ```

3. **Descargar los datos de secuenciación:**
   ```bash
   fastq-dump -X $N --split-files $SRR
   ```

4. **Alinear las lecturas con bwa y convertirlas a BAM:**
   ```bash
   bwa mem refs/$ACC.fa ${SRR}_1.fastq ${SRR}_2.fastq | samtools sort > $SRR.bwa.bam
   samtools index $SRR.bwa.bam
   ```

5. **Alinear las lecturas con bowtie2 (opcional):**
   ```bash
   bowtie2-build refs/$ACC.fa refs/$ACC.fa
   bowtie2 -x refs/$ACC.fa -1 ${SRR}_1.fastq -2 ${SRR}_2.fastq | samtools sort > $SRR.bowtie.bam
   samtools index $SRR.bowtie.bam
   ```

6. **Generar estadísticas del archivo BAM:**
   ```bash
   samtools flagstat $SRR.bwa.bam
   ```
---

### Ejemplo práctico: Análisis de un archivo SAM

Para comprender mejor los archivos SAM, se puede analizar un archivo utilizando herramientas como `samtools` y explorar las columnas clave.

1. **Ver el encabezado del archivo SAM:**
   ```bash
   samtools view -H $SRR.bwa.bam
   ```

2. **Extraer información de columnas específicas:**
   ```bash
   samtools view $SRR.bwa.bam | cut -f 1,3,4 | head -5
   ```
   Esto muestra los nombres de las lecturas, el nombre de la referencia y las posiciones de alineamiento.

3. **Explorar las etiquetas opcionales:**
   ```bash
   samtools view $SRR.bwa.bam | cut -f 12,13,14 | head -5
   ```
   Aquí se pueden analizar etiquetas como `NM` (distancia de edición) y `MD` (mismatches).

4. **Filtrar alineamientos por criterios específicos:**
   - Alineamientos en la hebra directa:
     ```bash
     samtools view -F 20 -b $SRR.bwa.bam > forward.bam
     ```
   - Alineamientos en la hebra inversa:
     ```bash
     samtools view -F 4 -f 16 -b $SRR.bwa.bam > reverse.bam
     ```

5. **Generar estadísticas detalladas:**
   ```bash
   samtools stats $SRR.bwa.bam > stats.txt
   ```
---

### Visualización con IGV

**Integrative Genomics Viewer (IGV)** es una herramienta ampliamente utilizada para visualizar alineamientos en formatos como BAM y CRAM.

1. **Descargar e instalar IGV:**
   - [Página oficial de IGV](https://software.broadinstitute.org/software/igv/)

2. **Preparar el archivo BAM/CRAM para IGV:**
   - Asegúrate de que el archivo está indexado:
     ```bash
     samtools index $SRR.bwa.bam
     ```

3. **Cargar los datos en IGV:**
   - Abre IGV y selecciona el genoma de referencia apropiado.
   - Ve a `File > Load from File` y selecciona el archivo BAM o CRAM.

4. **Navegar por los alineamientos:**
   - Usa la barra de búsqueda para ingresar una región específica (por ejemplo, `chr1:1000-2000`).
   - Explora las pistas de lectura para analizar variaciones, cobertura, y calidad de alineamiento.
