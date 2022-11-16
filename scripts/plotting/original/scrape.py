import os
import sys
import time
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains

lineage = sys.argv
del lineage[0]

if (not os.path.exists("./VariantSignatures")):
    os.mkdir("./VariantSignatures")

def get_mutations(lineage):

    options = Options()
    options.add_argument("--headless")

    driver = webdriver.Firefox(options=options)
    driver.set_window_size(3000,3000)

    link = "https://outbreak.info/compare-lineages?pango="+lineage+"&gene=ORF1a&gene=ORF1b&gene=S&gene=ORF3a&gene=E&gene=M&gene=ORF6&gene=ORF7a&gene=ORF7b&gene=ORF8&gene=N&gene=ORF10&threshold=50&nthresh=1&sub=false&dark=true"
    driver.get(link)

    time.sleep(60)

    if "." in lineage:
        match = lineage.replace(".","-")
        xpath = "//*[contains(@id,'"+match+"-')]"
    else:
        xpath = "//*[contains(@id,'"+lineage+"-')]"

    blocks = driver.find_elements(By.XPATH, xpath)

    file = "./VariantSignatures/"+lineage+".tsv"
    o = open(file, 'w')

    for element in blocks:
        hover = ActionChains(driver).move_to_element(element)
        hover.perform()
        percentage = driver.find_element(By.XPATH,"//div[@id='tooltip-prevalence'][contains(@style,'display: block')]//div[@id='prevalence']//div[@id='value']").text
        count = driver.find_element(By.XPATH,"//div[@id='tooltip-prevalence'][contains(@style,'display: block')]//div[@id='count']").text
        name = element.get_attribute("id")

        if "." in lineage:
            replacement = lineage.replace(".","-")+"-"
        else:
            replacement = lineage+"-"

        name = name.replace(replacement,"").replace("_",":",1).replace("_","-").upper()
        percentage = float(percentage.replace("%","")) / 100.0
        count = count.replace(')','').replace('(','')
        print(name,percentage,count)
        o.write("%s\t%.3f\t%s\n" % (name, percentage, count))
    o.close()
    driver.quit()

for element in lineage:
    print("Working on "+ element)
    get_mutations(element)

