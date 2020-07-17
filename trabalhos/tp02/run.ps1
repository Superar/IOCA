$files = Get-ChildItem -Name ".\data"

foreach ($file in $files) {
    Start-Job { python3 "D:\Documentos\Cursos\IOCA\trabalhos\tp02\solver.py" $args } -ArgumentList "D:\Documentos\Cursos\IOCA\trabalhos\tp02\data\$file"
}

Get-Job | Wait-Job
