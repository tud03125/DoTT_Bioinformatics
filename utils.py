def parse_conditions(bam_files, conditions):
    if conditions:
        cond_list = [cond.strip() for cond in conditions.split(",")]
        if len(cond_list) != len(bam_files):
            raise ValueError("Number of conditions provided does not match number of BAM files")
        return cond_list
    else:
        if len(bam_files) % 2 != 0:
            raise ValueError("Number of BAM files must be even if conditions are not provided")
        half = len(bam_files) // 2
        return ["Control"] * half + ["Experimental"] * half
