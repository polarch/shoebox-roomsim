function write_material(name,tag,absorption)
%WRITE_MATERIAL Summary of this function goes here
%   Detailed explanation goes here

material = [];
try
    load material_db
    disp('Existing material database "material_db.mat" found')
    Nmat = length(material);
    
    material(Nmat+1).name = name;
    material(Nmat+1).tag = tag;
    material(Nmat+1).absorption = absorption;
    if length(absorption) < 6
        error('Absorption coefficients should be defined in octaves from 125Hz - 4kHz or higher')
    end

    save material_db material    
    
catch
    disp('Existing material database not found')
    disp('Creating material database "material_db.mat"')
    material.name = name;
    material.tag = tag;
    material.absorption = absorption;
    if length(absorption) < 6
        error('Absorption coefficients should be defined in octaves from 125Hz - 4kHz or higher')
    end

    save material_db material
end

end

