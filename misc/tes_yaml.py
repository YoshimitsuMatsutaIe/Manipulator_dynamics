import yaml

with open('./misc/test.yaml') as file:
    config = yaml.safe_load(file.read())

print(config)