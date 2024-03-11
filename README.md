# SIBER-RISK Strong Motion Database

Este repositorio contiene los códigos y configuración necesarios para levantar la aplicación de corrección de registros sísmicos chilenos desarrollada en el marco del proyecto Fondecyt #1170836

## Estructura

```bash
.
├── data
│   └── seismicDatabase
│       ├── npz
│       └── mat
├── local
│   └── siberrisk
└── src
    ├── assets
    ├── data
    │   └── ffm
    ├── lib
    │   └── pyrjmcmc
    │      └── models
    │         ├── multiple_partition
    │         └── single_partition
    ├── secrets
    ├── tmp
    ├── websites
    │   ├── StrongMotionDatabase
    │   └── StrongMotionDatabaseDownloadManager
    └── widgets
```

## Código fuente

Todo el código relacionado tanto a la detección de la onda P, la corrección por línea base y la página web se encuentra en el directorio ```src```. Para utilizar la aplicación de corrección se debe ejecutar el archivo ```src/main.py```, por ejemplo desde la ruta ```src``` ejecutar el siguiente comando:

```shell
python3 main.py
````

En tanto, para levantar la página web del servidor se debe levantar un servidor local de ```bokeh``` con el archivo ```src/websites/StrongMotionDatabase/home.py```, por ejemplo desde la ruta ```src/websites``` ejecutar el siguiente comando:

```shell
bokeh serve StrongMotionDatabase/home.py --allow-websocket-origin=* --prefix=StrongMotionDatabase
```

La ejecución del código crea carpetas que no se sincronizan con GitHub necesariamente.

## Ejecución en local

Para tener certeza de poder levantar el sistema correctamente se ha creado un contenedor de Docker que permite ejecutar tanto la aplicación de corrección como la página web de forma local. Para esto, desde la ```local/siberrisk``` ejecutar los siguientes comandos:

```shell
docker compose build
docker compose up -d
```

Los servicios que se levantan se detallan en ```local/siberrisk/docker-compose.yaml```.

## Base de datos

La base de datos, ya sea corregida o no, se almacena en la ruta ```data```. Por razones de tamaño, esta carpeta no se sincroniza con GitHub. Los registros corregidos se guardan en archivos ```.npz``` para ser leidos en Python utilizando la librería NumPy, y en archivos ```.mat``` para ser leidos a través del programa MATLAB.

## Equipo

- [Sebastián Castro](https://github.com/sebacastroh)
- [Roberto Benavente](https://github.com/robenavente)
- Jorge G.F. Crempien
- [Gabriel Candia](https://github.com/gacandia)
- Juan Carlos de la Llera

## Acceso

La base de datos puede ser accedida libremente a través del enlace [https://doi.org/10.7764/datasetUC/ING-UC.1170836_1](https://doi.org/10.7764/datasetUC/ING-UC.1170836_1), en tanto la metodología puede ser consultada en el artículo científico asociado, disponible en [https://doi.org/10.1785/0220200336](https://doi.org/10.1785/0220200336).

## Cita

En caso de utilizar la base de datos, por favor utilice la siguiente cita

    Sebastián Castro, Roberto Benavente, Jorge G. F. Crempien, Gabriel Candia, Juan Carlos de la Llera; A Consistently Processed Strong‐Motion Database for Chilean Earthquakes. Seismological Research Letters 2022;; 93 (5): 2700–2718. doi: https://doi.org/10.1785/0220200336

y añada las siguientes líneas en la sección de agradecimientos

    The strong motion database was provided by the SIBER-RISK project: Simulation Based Earthquake Risk and Resilience of Interdependent Systems and Networks ANID/FONDECYT/1170836. doi: https://doi.org/10.7764/datasetUC/ING-UC.1170836_1
