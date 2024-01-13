#!/usr/bin/python
                    
"""                 
  Utility for SMATool -- Automated toolkit for computing zero and finite-temperature strength of materials
                    
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.
                
  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
                    
"""  

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options as ChromeOptions
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import os
import time
import gzip

# Configure the download directory
download_dir = 'qe_potentials'
if not os.path.exists(download_dir):
    os.makedirs(download_dir)

# Function to configure Chrome for automatic downloads
def configure_download_options(download_path):
    chrome_options = ChromeOptions()
    prefs = {"download.default_directory": download_path,
             "download.prompt_for_download": False,
             "download.directory_upgrade": True,
             "safebrowsing.enabled": True}
    chrome_options.add_experimental_option("prefs", prefs)
    return chrome_options

# Function to read options from pseudo.info
def read_pseudo_info():
    options = {
        "elements": None,
        "type": "NC SR (ONCVPSP v0.5)",
        "xc": "PBE",
        "accuracy": "standard",
        "format": "upf"
    }

    if os.path.exists("pseudo.info"):
        with open("pseudo.info", "r") as file:
            for line in file:
                if line.strip().lower().startswith('# list of elements to download'):
                    elements_line = next(file, '').strip()
                    options["elements"] = [elem.strip() for elem in elements_line.split(',')]
                elif "NC potential type" in line:
                    options["type"] = next(file, '').strip()
                elif "Select the DFT potential" in line:
                    options["xc"] = next(file, '').strip()
                elif "Accuracy of the DFT potential" in line:
                    options["accuracy"] = next(file, '').strip()
                elif "Format of the pseudopential" in line:
                    options["format"] = next(file, '').strip()
        print(f"The following elements and key info will be downloaded:\n{options}")
    else:
        print("pseudo.info not found. Using default options.")

    return options


pseudo_options = read_pseudo_info()
#print(pseudo_options["elements"])  

chrome_options = configure_download_options(os.path.abspath(download_dir))
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=chrome_options)

try:
    # Open the website
    driver.get("http://www.pseudo-dojo.org")
    time.sleep(5)  # Adjust based on page load time

    # Set the dropdown values
    Select(driver.find_element(By.ID, "TYP")).select_by_visible_text(pseudo_options["type"])
    Select(driver.find_element(By.ID, "XCF")).select_by_visible_text(pseudo_options["xc"])
    Select(driver.find_element(By.ID, "ACC")).select_by_visible_text(pseudo_options["accuracy"])
    Select(driver.find_element(By.ID, "FMT")).select_by_visible_text(pseudo_options["format"])
    time.sleep(5)  # Wait for the page to update

    # Identify clickable elements and initiate downloads
    elements = driver.find_elements(By.CSS_SELECTOR, '.plugin .element')
    for element in elements:
        element_name = element.text.strip()
        if pseudo_options["elements"] is None or element_name in pseudo_options["elements"]:
            element.click()
            time.sleep(1)  # Time for download to start

    # Wait for all downloads to complete
    time.sleep(30)  # Adjust based on the number of files and your internet speed

    # Rename downloaded files (if necessary)
    for filename in os.listdir(download_dir):
        if filename.endswith(".upf.gz"):
            # Extract element name from the filename
            element_name = filename.split('.')[0]
            new_filename = f"{element_name}_pdojo.upf"

            # Path for the compressed file and the new file
            gz_filepath = os.path.join(download_dir, filename)
            new_filepath = os.path.join(download_dir, new_filename)

            # Extract the .gz file
            with gzip.open(gz_filepath, 'rb') as f_in:
                with open(new_filepath, 'wb') as f_out:
                    f_out.write(f_in.read())

            # Optionally, remove the original .gz file
            os.remove(gz_filepath)

except Exception as e:
    print(f"An error occurred: {e}")

finally:
    # Close the browser
    driver.quit()



def main():
    pseudo_options = read_pseudo_info()

if __name__ == "__main__":
    main()
