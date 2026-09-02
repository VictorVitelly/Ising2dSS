#!/bin/bash
#SBATCH --job-name=Ising2dSSa                    # Nombre del trabajo
#SBATCH --output=logs/Ising2dSSa.log                   # Archivo de registro de salida
#SBATCH --error=logs/Ising2dSSa.err                    # Archivo de registro de errores
#SBATCH --partition=QuantPhysMC         # Nombre de la partición o cola de trabajos
#SBATCH --nodes=1                               # Número de nodos a utilizar (puedes cambiarlo)
#SBATCH --ntasks-per-node=1             # Número de tareas por nodo (1 para ejecución serial)
#SBATCH --cpus-per-task=1               # Número de CPUs por tarea (puedes cambiarlo)
#SBATCH --mem=4G                        # Memoria RAM necesaria (puedes cambiarlo)
#SBATCH --time=1-00

# Comando para ejecutar tu programa

cd ~/Ising2dSSa

SRC=src
BIN=bin
DAT=data
TARGET=main

(echo 1) |$BIN/$TARGET
