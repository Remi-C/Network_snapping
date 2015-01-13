-----------------------------------------------------------
--
--Rémi-C , Thales IGN
--01/2015
--
--This script converts open data paris trottoir to input for snapping
-----------------------------------------------------------

SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology, odparis_corrected, public; 


DROP TABLE IF EXISTS trottoir_cut_into_pieces ;  
CREATE TABLE trottoir_cut_into_pieces AS 
WITH trottoir AS ( --filtering to keep only ground cornerstones
	SELECT gid, geom 
	FROM trottoir
	WHERE niveau = 'Surface'
		AND libelle = 'Bordure'  
		AND ST_Area( Box2D(geom) ) >20
		AND ST_Length( geom) >14
)
SELECT row_number() over() as qgis_id, 1.0 as weight, 1.0 as confidence, gid, dmp.geom  as geom
FROM trottoir,ST_DumpPoints(ST_Segmentize(ST_SimplifyPreserveTopology(ST_Segmentize( geom, 4.0),0.5), 4.0)) AS dmp

CREATE INDEX ON trottoir_cut_into_pieces USING GIST (geom) ;
CREATE INDEX ON  trottoir_cut_into_pieces USING GIST (ST_Transform( geom,932011)) ; 
CREATE INDEX ON trottoir_cut_into_pieces (gid); 
--qgis_id,geom, weight ;

--cutting the result into isolated points, 
SELECT *
FROM weighted_sidewalk_observation_point
LIMIT 1