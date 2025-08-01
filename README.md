# Culturally Relevant Project-Based Learning and School Racial Climate

R code and analysis for a doctoral dissertation examining the impact of culturally relevant project-based learning on middle school students' perceptions of school racial climate.

## Study Overview

This quasi-experimental study compared three groups of middle school students:
- **Culturally Relevant PjBL Group**: Received project-based playwriting instruction with culturally relevant pedagogy
- **Standard PjBL Group**: Received standard project-based playwriting instruction  
- **Control Group**: Received traditional writing instruction

**Research Question**: Is there a difference among culturally relevant project-based writing instruction, standard project-based writing instruction, and traditional writing instruction in middle-school students' perceptions of school racial climate when controlling for prior English language arts achievement?

## Methods

- **Design**: Quasi-experimental posttest-only nonequivalent control group design
- **Participants**: 83 middle school students (grades 6 and 8) from urban Pennsylvania schools
- **Analysis**: One-way MANCOVA with planned contrasts
- **Instruments**: 
  - [School Climate for Diversity - Secondary Scale](https://byrdlab.wordpress.ncsu.edu/school-racial-climate/) (SCD-S; [Byrd](#references), 2017)
  - [Pennsylvania System of School Assessment](https://www.pa.gov/agencies/education/programs-and-services/instruction/elementary-and-secondary-education/assessment-and-accountability/pennsylvania-system-of-school-assessment-pssa) (PSSA) ELA test ([Pennsylvania Department of Education](#references))

## Repository Contents

- `analysis.md` - **Start Here** GitHub-rendered analysis with full output and explanations
- `analysis.Rmd` - Annotated R Markdown file 
- `analysis.R` - Raw R code without annotations
- Data files are not included to maintain participant privacy

## Key Findings

The omnibus MANCOVA found no statistically significant differences between groups (p > .05). However, planned contrasts revealed that combined PjBL groups showed significantly higher scores than the control group on:
- Promotion of cultural competence (p = .015, d = 0.65)
- Critical consciousness socialization (p = .022, d = 0.62)

## Software Requirements

- R (version 4.5.1 or later)
- Required packages listed in the R Markdown file

## Acknowledgments

Some portions of the R code were developed with assistance from AI tools including Claude.ai and Gemini for coding assistance and troubleshooting.

## Citation

If you use this code in your research, please cite:

Aungst, G. W. (2025). *Culturally relevant project-based learning analysis code* [Computer software]. GitHub. https://github.com/geraldaungst/pjbl-dissertation

For the full dissertation, cite:

Aungst, G. W. (2025). *Culturally relevant project-based learning: A quasi-experimental study of the effect of playwriting instruction on student perceptions of school racial climate* [Doctoral dissertation, Liberty University].

## References

Byrd, C. M. (2017). The complexity of school racial climate: Reliability and validity of a new measure for secondary students. *British Journal of Educational Psychology*, *87*(4), 700-721. https://doi.org/10.1111/bjep.12179

Pennsylvania Department of Education. (n.d.). Pennsylvania System of School Assessment (PSSA). Retrieved from https://www.education.pa.gov/K-12/Assessment%20and%20Accountability/PSSA/Pages/default.aspx

## Contact

Gerald W. Aungst  
Liberty University  
[gwaungst@liberty.edu](mailto:gwaungst@liberty.edu)

## License

This project is licensed under the MIT License - see the LICENSE file for details.
