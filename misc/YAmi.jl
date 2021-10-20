import YAML
data = YAML.load_file("./config./test.yaml")
println(data["rmp_param"])