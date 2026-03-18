# Load the molecule .pdb file
mol new [structure file]
mol addrep 0
# Atom settings
mol modcolor 0 0 Element
mol modmaterial 0 0 Glossy
mol modstyle 0 0 CPK 0.8 0.0 25.0 1.0
# Bond settings
mol modcolor 1 0 Element
mol modmaterial 1 0 Glossy
mol modstyle 1 0 DynamicBonds 1.560000 0.03500000 300.000000

# Load CPs from CPs.pdb
mol new [CPs file]
# Loading creates Rep 0 (Lines all). We MODIFY it.

# --- CP VISUALIZATION (BACKGROUND - HIDDEN) ---

# Nuclear (3,-3) - C - Rep 0
mol modcolor 0 1 ColorID 6  ;# Silver
mol modselect 0 1 element C
mol modstyle 0 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 0 1 Glossy
mol showrep 1 0 0 ;# Hidden

# Bond (3,-1) - N - Rep 1
mol addrep 1
mol modcolor 1 1 ColorID 12 ;# Lime
mol modselect 1 1 element N
mol modstyle 1 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 1 1 Glossy
mol showrep 1 1 0 ;# Hidden

# Ring (3,+1) - O - Rep 2
mol addrep 1
mol modcolor 2 1 ColorID 11 ;# Purple
mol modselect 2 1 element O
mol modstyle 2 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 2 1 Glossy
mol showrep 1 2 0 ;# Hidden

# Cage (3,+3) - F - Rep 3
mol addrep 1
mol modcolor 3 1 ColorID 22 ;# Yellow/Green
mol modselect 3 1 element F  
mol modstyle 3 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 3 1 Glossy
mol showrep 1 3 0 ;# Hidden


# --- INTEREST CP VISUALIZATION (VISIBLE) ---

# INTEREST - Nuclear (C) - Rep 4
mol addrep 1
mol modcolor 4 1 ColorID 6 ;# Silver
mol modselect 4 1 "index [interest_cps] and element C"
mol modstyle 4 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 4 1 Glossy
mol showrep 1 4 1 

# INTEREST - Bond (N) - Rep 5
mol addrep 1
mol modcolor 5 1 ColorID 12 ;# Lime
mol modselect 5 1 "index [interest_cps] and element N"
mol modstyle 5 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 5 1 Glossy
mol showrep 1 5 1 

# INTEREST - Ring (O) - Rep 6
mol addrep 1
mol modcolor 6 1 ColorID 11 ;# Purple
mol modselect 6 1 "index [interest_cps] and element O"
mol modstyle 6 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 6 1 Glossy
mol showrep 1 6 1 

# INTEREST - Cage (F) - Rep 7
mol addrep 1
mol modcolor 7 1 ColorID 22 ;# Yellow/Green
mol modselect 7 1 "index [interest_cps] and element F"
mol modstyle 7 1 CPK 0.3 0.0 50.0 1.0
mol modmaterial 7 1 Glossy
mol showrep 1 7 1 


# --- PATH VISUALIZATION ---
# Load paths.pdb (Molecule 2)
mol new [paths file] filebonds off autobonds off
# Rep 0 (Lines all) is created by default.

# Rep 0: General Paths (Background) -> Hide
mol modstyle 0 2 Lines 0.1
mol modcolor 0 2 ColorID 11
mol modmaterial 0 2 Transparent
mol modselect 0 2 all
mol showrep 2 0 0 ;# HIDDEN

# Rep 1: Interest Paths (Highlighted) -> Visible
mol addrep 2
mol modstyle 1 2 Lines 2.0
mol modcolor 1 2 ColorID 11
mol modmaterial 1 2 Opaque
mol modselect 1 2 "[selection_paths]"
mol showrep 2 1 [show_paths_rep]

# Global Visual Settings
color Element C silver
color Element N blue
color Element H white
color Element O red
color Element F yellow3
color Element Cl green

display resetview
color Display Background white
