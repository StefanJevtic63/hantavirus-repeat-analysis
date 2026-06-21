@echo off

:: Force the terminal to focus on the directory where this .bat file lives
cd /d "%~dp0"

:: Check if the executable exists in this directory
if not exist StatRepeatsNoDB.exe (
    echo Error: StatRepeatsNoDB.exe not found in the current directory!
    echo Current directory is: %cd%
    pause
    exit /b 1
)

:: --------------------------------------------------
:: Cleanup Phase 
:: --------------------------------------------------
echo Cleaning up and recreating 'data\output' directory...
if exist data\output rmdir /s /q data\output
mkdir data\output

echo Preparing input files...
python preprocess\clean_fasta.py


:: --------------------------------------------------
:: Analysis Phase
:: --------------------------------------------------
echo.
echo --------------------------------------------------
echo Starting analysis...
echo --------------------------------------------------

:: Process: all_sequences_clean.fasta
echo Processing: all_sequences_clean dn...
StatRepeatsNoDB.exe data\input\all_sequences_clean.fasta 5 -rna -dn -output data\output\all_dn.txt

echo.
echo Processing: all_sequences_clean dc...
StatRepeatsNoDB.exe data\input\all_sequences_clean.fasta 5 -rna -dc -output data\output\all_dc.txt

echo.
echo Processing: all_sequences_clean in...
StatRepeatsNoDB.exe data\input\all_sequences_clean.fasta 5 -rna -in -output data\output\all_in.txt

echo.
echo Processing: all_sequences_clean ic...
StatRepeatsNoDB.exe data\input\all_sequences_clean.fasta 5 -rna -ic -output data\output\all_ic.txt

:: Process: complete_sequences_clean.fasta
echo.
echo Processing: complete_sequences_clean dn...
StatRepeatsNoDB.exe data\input\complete_sequences_clean.fasta 5 -rna -dn -output data\output\complete_dn.txt

echo.
echo Processing: complete_sequences_clean dc (DNA mode fallback)...
StatRepeatsNoDB.exe data\input\complete_sequences_clean.fasta 5 -dc -output data\output\complete_dc.txt

echo.
echo Processing: complete_sequences_clean in...
StatRepeatsNoDB.exe data\input\complete_sequences_clean.fasta 5 -rna -in -output data\output\complete_in.txt

echo.
echo Processing: complete_sequences_clean ic...
StatRepeatsNoDB.exe data\input\complete_sequences_clean.fasta 5 -rna -ic -output data\output\complete_ic.txt


:: --------------------------------------------------
:: CSV Creation Phase
:: --------------------------------------------------
echo.
echo Creating the combined CSV file from the analysis outputs...
python preprocess\create_csv.py
