pyinstaller --onefile --name msmutect  --hidden-import=tempfile --hidden-import=json --hidden-import=ctypes ../src/Entry/main.py