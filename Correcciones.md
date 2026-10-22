# Correcciones

Revisión del proyecto `dna-compression-llm`. El código funcionaba en lo básico,
pero tenía fallos que las pruebas existentes no detectaban. Todos se
reprodujeron antes de corregirlos y ahora tienen una prueba que los cubre.

## Fallos corregidos

1. **La restauración no era idéntica al original.** `compress` eliminaba las
   líneas vacías y los espacios al inicio y al final de cada línea, aunque el
   proyecto promete una copia byte a byte. Ahora solo se separa el salto de
   línea y se conserva todo lo demás.

2. **Algunos archivos comprimidos no se podían recuperar.** Si el nombre de un
   token contenía letras que coincidían con un patrón posterior, ese patrón
   reemplazaba parte del token (por ejemplo `<<#>>`) y la descompresión
   fallaba. Ahora los reemplazos solo se aplican al ADN, nunca a tokens ya
   insertados. Además:
   - se valida el JSON de patrones (tokens y secuencias no vacíos y sin `<` ni `>`);
   - `compress` rechaza con un error claro las líneas que contienen `<` o `>`.

3. **El ahorro estimado era incorrecto.** `detect` calculaba el ahorro como
   `largo - 1`, pero cada reemplazo escribe `<X>`, que ocupa 3 bytes. Por eso
   conservaba patrones que hacían crecer el archivo, como `GG`. Ahora `detect` y
   `analyze` usan la misma fórmula: `count × (largo - 3)`.

4. **Las ocurrencias solapadas se contaban de más.** `ATAT` en `ATATATATAT`
   contaba 4, pero el compresor solo puede hacer 2 reemplazos. Ahora se cuentan
   solo las ocurrencias que no se solapan.

5. **`benchmark` fallaba con un archivo vacío** (división entre cero).

## Mejoras

- **API keys:** se leen de `OPENAI_API_KEY`, `DEEPSEEK_API_KEY` o
  `GEMINI_API_KEY`, y `--key` sigue funcionando. Antes, sin `--key`, se enviaba
  `"local"` como clave y cada lote fallaba sin avisar. Ahora el programa se
  detiene con un mensaje claro.
- **Gemini** recibe las instrucciones de sistema y `temperature=0`, igual que
  ChatGPT y DeepSeek.
- **`--model` con modelos locales** se pasa directamente al cargar el modelo,
  sin modificar la configuración global.
- **Errores legibles:** si falta un archivo o el JSON es inválido, se muestra
  un mensaje en vez de una traza de Python.
- **`extract`** rechaza archivos que no son FASTQ o que están truncados.
- **Log:** se muestra la ruta real del log y se registra también la carga del modelo.
- **Validación de argumentos:** `--batch_size`, `--threads`, `--step`,
  `--max-overhead`, `--top` y `--n` deben ser mayores que cero.
- **Documentación:** se actualizaron `README.md` y `docs/architecture.md`.

## Compatibilidad

- Los JSON de patrones y los archivos `.compress` de la versión anterior siguen
  funcionando.
- Con el mismo diccionario, la versión anterior y la nueva producen el mismo
  archivo comprimido.
- `detect` puede generar un diccionario algo distinto que antes, porque ahora
  descarta patrones que no ahorran espacio.

## Verificación

- 45 pruebas pasan (antes 31) con Python 3.10.
- Flujo completo con `data/input/ERR15993673_5000.seq.txt`: el archivo
  restaurado tiene el mismo SHA-256 que el original.
- `uv lock --check` y `uv build` funcionan sin cambios en las dependencias.

## Pendiente (no se modificó)

- El proveedor `hyenadna` calcula embeddings y luego no los usa. Los patrones
  salen de una búsqueda de repeticiones. `dnabert2` usa sobre todo el tokenizador.
- Varias funciones están duplicadas entre `overhead.py` y `detection/pipeline.py`.
- Con los datos de ejemplo la compresión es de alrededor del 2 %, mientras que
  gzip logra un 95 % y lzma un 97 %.
- No se probó con APIs ni modelos reales; esa parte está cubierta con
  proveedores simulados.

## Archivos modificados

- `README.md`
- `docs/architecture.md`
- `dna_compression/benchmark.py`
- `dna_compression/cli.py`
- `dna_compression/compression.py`
- `dna_compression/config.py`
- `dna_compression/fastq.py`
- `dna_compression/overhead.py`
- `dna_compression/pattern_store.py`
- `dna_compression/detection/pipeline.py`
- `dna_compression/detection/providers.py`
- `tests/test_cli.py`
- `tests/test_compression.py`
- `tests/test_detection_pipeline.py`
- `tests/test_fastq.py`
- `tests/test_providers.py`
- `tests/test_utilities.py`
