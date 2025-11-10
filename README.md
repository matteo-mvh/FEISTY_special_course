# FEISTY Special Course

**FEISTY_special_course** explores marine ecosystem dynamics and the **impact of fishing on the biological carbon pump** using the **FEISTY model**.  
The repository is divided into two main components: **Local** and **Global**, which together bridge small-scale ecosystem understanding and global-scale carbon impact assessment.

---

## Repository Structure

### 🐟 Local
The **Local** section focuses on **specific marine regions** to study ecosystem processes and fishing impacts under controlled conditions.  
It includes three ecosystems differing in **depth** and **productivity**:

- **Shelf Sea (75 m)** – highly productive, dominated by demersal and forage fish.  
- **Slope (1500 m)** – intermediate depth with mixed trophic interactions.  
- **Open Ocean (3000 m)** – low productivity, limited biomass.

**Objectives**
- Explore equilibrium **biomass distribution** and **food-web structure** at local scales.  
- Test the effects of **size-selective and targeted fishing** (demersal, forage, large pelagic).  
- Estimate **carbon injection potential** for different fishing scenarios.  

**Outputs**
- Equilibrium **biomass** per size class and depth.  
- **Carbon injection** by fishing type and its **difference to the baseline**.  
- Visuals showing **biomass shifts** and **carbon export efficiency**.

---

### 🌍 Global
The **Global** section runs FEISTY at a **planetary scale**, integrating environmental forcing data and a **Global Transport Matrix** to link biological production with carbon storage.

**Objectives**
- Quantify **total global biomass** and **carbon fluxes**.  
- Assess the **global effects of fishing** on carbon sequestration.  
- Combine biological and physical models to evaluate **long-term carbon storage**.

**Outputs**
- **Global biomass maps** and **fishing response** visualizations.  
- **Carbon injection and sequestration estimates** per fishing scenario.  
- Aggregated plots showing **differences from unfished baseline** conditions.

---

## Key Features

- Modular and scalable **FEISTY simulations** (Local → Global).  
- Integration of **size-structured fishing** and **ecosystem depth variation**.  
- Calculation of **biomass**, **carbon injection**, and **carbon sequestration**.  
- Outputs formatted for **visualization** and **further data analysis**.  

---

## Usage

1. Clone the repository:  
   ```bash
   git clone https://github.com/matteo-mvh/FEISTY_special_course.git
