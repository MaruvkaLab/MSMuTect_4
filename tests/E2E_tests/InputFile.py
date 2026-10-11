from dataclasses import dataclass

@dataclass
class InputFile:
    local_file_directory: str
    docker_file_directory: str
    filename: str
    flag: str

