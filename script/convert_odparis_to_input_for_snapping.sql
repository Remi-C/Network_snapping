-----------------------------------------------------------
--
--Rémi-C , Thales IGN
--01/2015
--
--This script converts open data paris trottoir to input for snapping
-----------------------------------------------------------

SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology, odparis_corrected, public; 

/*
 DROP TABLE IF EXISTS user_defined_sidewalk_border ;  
CREATE TABLE user_defined_sidewalk_border ( 
gid serial primary key
, weight float DEFAULT 500
, geom geometry(Multipoint, 932011) 
);
CREATE INDEX ON user_defined_sidewalk_border USING GISt(geom) ; 
CREATE INDEX ON user_defined_sidewalk_border  (weight); 


*/
DROP TABLE IF EXISTS trottoir_cut_into_pieces ;  
CREATE TABLE trottoir_cut_into_pieces AS 
WITH trottoir AS ( --filtering to keep only ground cornerstones
	SELECT t.gid, t.geom 
	FROM trottoir as t  , def_zone_export as dfz 
	WHERE niveau = 'Surface'
		AND libelle = 'Bordure'  
		AND ST_Area( Box2D(t.geom) ) >20
		AND ST_Length( t.geom) >14
		AND ST_INtersects( ST_Transform(t.geom,932011), ST_Transform(dfz.geom,932011))=TRUE 
)
SELECT row_number() over() as qgis_id, 1.0 as weight, 1.0 as confidence, gid, dmp.geom  as geom
FROM trottoir,ST_DumpPoints(ST_Segmentize(ST_SimplifyPreserveTopology(ST_Segmentize( geom, 4.0),0.5), 4.0)) AS dmp ; 

CREATE INDEX ON trottoir_cut_into_pieces USING GIST (geom) ;
CREATE INDEX ON  trottoir_cut_into_pieces USING GIST (ST_Transform( geom,932011)) ; 
CREATE INDEX ON trottoir_cut_into_pieces (gid); 
--qgis_id,geom, weight ;
 

 
DROP TABLE IF EXISTS obs_for_output_in_export_area; 
	CREATE TABLE obs_for_output_in_export_area AS  
	WITH edge_geom AS ( -- we reconstruct the edge geom to be able to assign observation to edges
		SELECT efo.*, ST_SetSRID(ST_MakeLine(ST_MakePoint(nfo1.X,nfo1.Y,nfo1.Z) ,ST_MakePoint(nfo2.X,nfo2.Y,nfo2.Z)  ),932011) as edge_geom
		FROM edges_for_output_in_export_area as efo
			LEFT JOIN nodes_for_output_in_export_area AS nfo1 ON (efo.start_node = nfo1.node_id)
			LEFT JOIN nodes_for_output_in_export_area AS nfo2 ON (efo.end_node = nfo2.node_id)
	)
	,map AS (--for each edge, we get observation closer than 5 meters
		SELECT DISTINCT ON (oia.qgis_id) ST_Transform(oia.geom,932011) as geom, oia.weight, eg.edge_id
		FROM trottoir_cut_into_pieces   AS oia  ,  edge_geom AS eg
		WHERE ST_DWithin( ST_Transform(oia.geom,932011), eg.edge_geom,4+width)=TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface,ST_Transform(oia.geom,932011),10)=TRUE
			)
		ORDER BY oia.qgis_id ASC, ST_Distance(ST_Transform(oia.geom,932011),eg.edge_geom ) ASC
			
	)
	,map2 AS (
		SELECT DISTINCT ON (oia.gid) dmp.geom, oia.weight, eg.edge_id
		FROM user_defined_sidewalk_border   AS oia  ,  edge_geom AS eg, ST_DumpPoints(geom) AS dmp
		WHERE ST_DWithin(  oia.geom , eg.edge_geom,4+width)=TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface, oia.geom ,10)=TRUE
			)
		ORDER BY oia.gid ASC, ST_Distance( oia.geom ,eg.edge_geom ) ASC
	)
	SELECT row_number() over() AS obs_id , edge_id
		, ST_X(geom) AS X
		, ST_Y(geom) AS Y
		, COALESCE(ST_Z(geom),0) AS Z 
		, 1::float AS confidence, COALESCE(weight,0) AS weight
	FROM (SELECT * FROM map UNION ALL SELECT * FROM map2) as umap  ; 

