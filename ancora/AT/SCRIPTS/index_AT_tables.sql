-- Index MAPPING table
ALTER TABLE MAPPING ADD INDEX (qName);
ALTER TABLE MAPPING ADD INDEX chr_start (tName, bin, tStart);
ALTER TABLE MAPPING ADD INDEX chr_end (tName, bin, tEnd);

-- Index HSP table
ALTER TABLE HSP ADD INDEX (mapping_id);
ALTER TABLE HSP ADD INDEX (tStart);
ALTER TABLE HSP ADD INDEX (tEnd);

--- Index mrna table
--- Note: before doing this you should remove entries from the mrna
--- table with no counterpart in the MAPPING table; this greatly reduces
--- the size of the mrna table
ALTER TABLE mrna ADD UNIQUE KEY (acc);
ALTER TABLE mrna ADD KEY (cds);
ALTER TABLE mrna ADD KEY (mrnaClone);
ALTER TABLE mrna ADD KEY (library);

--- Index auxilary tables
ALTER TABLE cds ADD PRIMARY KEY (id);
ALTER TABLE mrnaClone ADD PRIMARY KEY (id);
ALTER TABLE library ADD PRIMARY KEY (id);
ALTER TABLE refSeqStatus ADD PRIMARY KEY (mrnaAcc);

