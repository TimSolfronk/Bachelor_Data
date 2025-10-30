with open("dependencies.txt", "r") as deps:
    package_names = (deps.read()).split("^")
    install_command = ""
    for package in package_names:
        package_name_with_version = package.split(" ")[0]
        install_command += "^" + package_name_with_version + " "

    install_command = "spack install " + install_command[1:]
    print(install_command)

