from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
import time
import sys

# The URL of the page with the form
STRIKE_URL = "https://strike.scec.org/cvws/cgi-bin/cvws.cgi"
scenario_name = "tpv5"
reference_name = "barall"

PUB_AREA_CODE = "G0012"
SCENARIOS_CODE = "G1045"
SOLUTIONS_CODE = "G1047"
RAW_DATA_CODE = "G1059"


def get_all_possible_options(option_code):
    option_names = []
    curIndex = -1
    rest_of_page = driver.page_source
    nextIndex = rest_of_page[1:].find(option_code)
    while nextIndex != -1:
        curIndex = nextIndex+1
        rest_of_page = rest_of_page[curIndex:]
        until = rest_of_page.find('"')
        new_option = rest_of_page[len(option_code):until]
        #print(new_option)
        option_names.append(new_option)

        nextIndex = rest_of_page[1:].find(option_code)

    return option_names

#format should be: tracer_id 0 t x y z data
#input format of curline is: t data
def format_data_line(tracer_id, tracer_pos, cur_line_data: str):
    # TODO: Format properly
    line = str(tracer_id)+",0"
    data_entries = cur_line_data.strip().split(" ")
    
    line += "," + data_entries[0].strip()
    
    for coord in tracer_pos:
        line += "," + str(coord)
    for i in range(1,len(data_entries)):
        if data_entries[i] == "":
            continue
        line += "," + data_entries[i].strip()

    return line + "\n"

def get_single_coordinate(tracer_name):
    if tracer_name.startswith("-"):
        return - float(tracer_name[1:3] + "." + tracer_name[3])
    else:
        return float(tracer_name[:2] + "." + tracer_name[2])

# converts a tracer_name that is either structured as 'faultst-045dp000' or 'body030st-120dp075' to the corresponding tracer position
def get_tracer_pos(tracer_name:str):
    x = 0.0
    if tracer_name.startswith("fault"):
        x = 20.0
        tracer_name = tracer_name[5:]
    elif tracer_name.startswith("body"):
        tracer_name = tracer_name[4:]
        x = get_single_coordinate(tracer_name)
        tracer_name = tracer_name[3+(1 if x < 0.0 else 0):]
        x += 20.0
    
    #skip 'st'
    tracer_name = tracer_name[2:]
    z = get_single_coordinate(tracer_name)
    tracer_name = tracer_name[3+(1 if z < 0.0 else 0):]
    z += 20.0

    #skip 'dp'
    tracer_name = tracer_name[2:]
    y = get_single_coordinate(tracer_name)

    return (x,y,z)

opts = Options()
# opts.add_argument("--headless=new")   # uncomment to run headless (no visible browser)
driver = webdriver.Chrome(options=opts)  # ensure chromedriver is installed/available

driver.get(STRIKE_URL)  # page that contains the form
# find the submit input by name and click it
pub_area_button = driver.find_element(By.NAME, PUB_AREA_CODE)
pub_area_button.click()

# give the browser a moment to follow redirect / load
time.sleep(1)
try:
    scenario_button = driver.find_element(By.NAME,SCENARIOS_CODE + scenario_name)
    scenario_button.click()
    time.sleep(1)
except:
    print("Couldn't find scenario '" + scenario_name + "'.")
    print("Available scenarios are: " + str(get_all_possible_options(SCENARIOS_CODE)))
    print("Exiting now.")
    driver.quit()
    sys.exit()

try:
    ref_sol_button = driver.find_element(By.NAME,SOLUTIONS_CODE + reference_name)
    ref_sol_button.click()
    time.sleep(1)
except:
    print("Couldn't find reference solution '" + reference_name + "' for scenario '" + scenario_name + "'.")
    print("Available reference solutions are: " + str(get_all_possible_options(SOLUTIONS_CODE)))
    print("Exiting now.")
    driver.quit()
    sys.exit()

tracer_names = get_all_possible_options(RAW_DATA_CODE)

with open(scenario_name.upper() + "_ref_" + reference_name + ".csv", "w") as output_file:
    
    output_file.write("number(0), number(1), t, x(0), x(1), x(2), data \n")
    for i in range(0,len(tracer_names)):
        raw_data_button = driver.find_element(By.NAME,RAW_DATA_CODE + tracer_names[i])
        raw_data_button.click()
        time.sleep(1)
        raw_data = driver.find_element(By.TAG_NAME, "pre").text
        raw_data_lines = raw_data[raw_data.find("\nt ")+1:].split("\n")
        raw_data_lable = raw_data_lines[0]
        tracer_pos = get_tracer_pos(tracer_names[i])

        for j in range(1,len(raw_data_lines)):
            if raw_data_lines[j].startswith("#"):
                continue
            output_file.write(format_data_line(i,tracer_pos,raw_data_lines[j]))

        driver.back()
        time.sleep(0.5)

#print(tracer_names)

#print(driver.page_source[:])
driver.quit()
