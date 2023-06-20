import os


class IKineticReader:
    def read(self, file_path: str) -> dict[str, tuple[float, float]]:
        raise NotImplementedError("Subclass must implement read method")


class CSVReader(IKineticReader):
    def read(self, file_path: str) -> dict[str, tuple[float, float]]:
        print("Reading CSV file")


class TXTReader(IKineticReader):
    def read(self, file_path: str) -> dict[str, tuple[float, float]]:
        print("Reading TXT file")


class TSVReader(IKineticReader):
    def read(self, file_path: str) -> dict[str, tuple[float, float]]:
        print("Reading TSV file")


class FileReaderFactory:
    def get_file_reader(self, file_path):
        _, file_extension = os.path.splitext(file_path)
        if file_extension == ".csv":
            return CSVReader()
        elif file_extension == ".txt":
            return TXTReader()
        elif file_extension == ".tsv":
            return TSVReader()
        else:
            raise ValueError("Invalid file extension")
