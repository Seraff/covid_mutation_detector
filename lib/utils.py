def create_dir_if_not_exists(path):
  Path(path).mkdir(parents=True, exist_ok=True)
