-----------------------------------------------------------
--
--Rémi-C , Thales IGN
--09/2014
--
--This script converts detected objects to input for snapping : 
--#object_id;class_id;class_name;edge_id;geom;confidence;
--
----- 
--warning : uses SFCGAL
-----------------------------------------------------------


CREATE SCHEMA IF NOT EXISTS objects_for_snapping; 

SET search_path TO objects_for_snapping,network_for_snapping, bdtopo_topological, bdtopo, topology, public; 


	--def_zone_export  
 
  /* --this part are only to be done one time
	--import table with objects :
	DROP TABLE IF EXISTS ashaped_objects ; 
	CREATE TABLE ashaped_objects AS
		SELECT *
		FROM public.dblink(
		'hostaddr=127.0.0.1 port=5433 dbname=TerraMobilita user=postgres password=postgres'::text
		, 'SELECT gid, classification_id,classification,z_range ,geom,import_error,weighted_confidence,should_be_used 
		FROM cmm_alpha_shaped_objects.ashaped_objects'::text
		,TRUE) AS
		f(
			gid INTEGER,
			  classification_id integer,
			  classification text,  
			  z_range numrange,
			  geom geometry(MultiPolygon,932011), 
			  import_error boolean,
			  weighted_confidence double precision,
			  should_be_used double precision
		) ;
		CREATE INDEX ON ashaped_objects USING GIST(geom) ; 
		CREATE INDEX ON ashaped_objects  (classification_id) ; 
		CREATE INDEX ON ashaped_objects  (classification) ; 
		CREATE INDEX ON ashaped_objects USING GIST (z_range) ; 
		CREATE INDEX ON ashaped_objects  (import_error) ; 
		CREATE INDEX ON ashaped_objects  (weighted_confidence) ; 
		CREATE INDEX ON ashaped_objects  (should_be_used) ; 
		CREATE INDEX ON ashaped_objects  (gid) ; 
		--ANALYZE ashaped_objects;
		--REFRESH MATERIALIZED VIEW ashaped_objects ; 
 
	--mapping the object to the network 
	

	--REFRESH MATERIALIZED VIEW ashaped_objects ;
*/
/*
	DROP TABLE IF EXISTS obj_for_output_in_export_area; 
	CREATE TABLE obj_for_output_in_export_area AS 
	WITH obj_in_area AS ( --these are the observations that are inside the defined area
-- 			SELECT  ao.gid ,  classification_id ,  classification ,  ST_Envelope(ao.geom) AS geom , weighted_confidence ,   should_be_used 
-- 			FROM ashaped_objects AS  ao, def_zone_export as dfz
-- 			WHERE ST_WITHIN(ao.geom,  ST_Transform(dfz.geom ,932011) )= TRUE 
-- 				AND (upper(z_range)-lower(z_range)) BETWEEN 0.2 AND 4 
-- 				AND classification = ANY (ARRAY[
-- 					'bollard'
-- 					,'other object'
-- 					,'pedestrian'
-- 					,'parked bicycle'
-- 					,'traffic sign'
-- 					--,'sidewalk'
-- 					,'walking pedestrian'
-- 					--,'car'
-- 					--,'other'
-- 					--,'tree'
-- 					,'trash can'
-- 					,'holding pedestrian'
-- 					,'scooter without driver'
-- 					--,'no_classification'
-- 					--,'building'
-- 					--,'road'
-- 					,'meter'
-- 					--,'ground'
-- 					,'still pedestrian'
-- 					,'punctual object'
-- 					--,'curb'
-- 					,'potted plant' 
-- 					])
-- 			UNION ALL
			SELECT ao.gid ,  classification_id ,  classification 
					,  ST_SnapToGrid(ST_SimplifyPreserveTopology(ST_Buffer(ST_Buffer(ao.geom ,2,'quad_segs=4'),-2, 'quad_segs=4'),0.1),0.001) AS geom
					, weighted_confidence ,   should_be_used 
			FROM ashaped_objects AS  ao, def_zone_export as dfz
			WHERE ST_WITHIN(ao.geom,  ST_Transform(dfz.geom ,932011) )= TRUE 
				AND classification = 'car'
				AND (upper(z_range)-lower(z_range)) BETWEEN 0.2 AND 4 
				AND ST_Area(ao.geom) BETWEEN 3.0 AND 12.0
		)
	,edge_geom AS ( -- we reconstruct the edge geom to be able to assign objects to edges
		SELECT efo.*, ST_SetSRID(ST_MakeLine(ST_MakePoint(nfo1.X,nfo1.Y,nfo1.Z) ,ST_MakePoint(nfo2.X,nfo2.Y,nfo2.Z)  ),932011) as edge_geom
		FROM edges_for_output_in_export_area as efo
			LEFT JOIN nodes_for_output_in_export_area AS nfo1 ON (efo.start_node = nfo1.node_id)
			LEFT JOIN nodes_for_output_in_export_area AS nfo2 ON (efo.end_node = nfo2.node_id)
	)
	,map AS (--for each edge, we get objects closer than 5 meters
		SELECT 
			DISTINCT ON (oia.gid) --note : enable = object mapped to at most 1 edge
			oia.*, eg.edge_id
			,ST_ShortestLine(ST_Transform( oia.geom ,932011),  eg.edge_geom) as sline_to_edge
		FROM obj_in_area  AS oia ,edge_geom AS eg
		WHERE ST_DWithin( ST_Transform(oia.geom,932011), eg.edge_geom,4+width)=TRUE
			--order by usefull if using distinct
				AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface,ST_Transform(oia.geom,932011),10)=TRUE
			)
		ORDER BY oia.gid ASC, abs(ST_Distance(ST_Transform(oia.geom,932011),ST_Buffer(eg.edge_geom,width/2,'endcap=flat' ))) ASC 
	)
	--object_id;class_id;class_name;edge_id;geom;confidence;
	SELECT  
		row_number() over() AS qgis_id
		,gid AS object_id 
		,classification_id AS class_id
		, classification AS class_name
		,   edge_id 
		, ST_Force2D(geom)  AS geom 
		, ST_Force2D(sline_to_edge) AS sline_to_edge
		,weighted_confidence AS confidence 
	FROM map  ; 
*/

 
	

--trying to accelerate the querry : it doesn't scale ! 

	DROP TABLE IF EXISTS obj_for_output_in_export_area; 
	CREATE TABLE obj_for_output_in_export_area AS   
	WITH crude_mapping AS ( --it is very counter intuitive, but it isway faster to force the data base to use the geometric filtering on inner product first, then use other
		SELECT   oia.gid
			,oia.classification
			,oia.z_range
			,oia.classification_id
			,oia.weighted_confidence
			,oia.geom
			, eg.edge_id
			, eg.width 
			,eg.geom AS edge_geom
		FROM ashaped_objects   AS oia, edges_for_output AS eg , def_zone_export as dfz
		WHERE  
		 ST_DWithin(  oia.geom , eg.geom,4+width) = true
	)
	,map AS (
		SELECT DISTINCT ON (oia.gid) oia.*
			,ST_ShortestLine( oia.geom  , edge_geom) as sline_to_edge
		FROM crude_mapping as oia
			INNER JOIN edges_for_output_in_export_area USING(edge_id) --enforce the fact that all edge are going to be exported
			, def_zone_export as dfz 
		WHERE classification = 'car'
			AND (upper(z_range)-lower(z_range)) BETWEEN 0.2 AND 4 
			AND ST_Area( oia.geom) BETWEEN 3.0 AND 12.0
			AND  ST_WITHIN(oia.geom,  ST_Transform(dfz.geom ,932011) )= TRUE  --keeping only object that are in the export are 
			AND NOT EXISTS ( --checking that object are not too close to intersection, or we don't consider theim
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface, oia.geom ,10)=TRUE
			)    
		--order by usefull if using distinct
		ORDER BY  oia.gid ASC, abs(ST_Distance(  oia.geom ,ST_Buffer(edge_geom,oia.width/2,'endcap=flat' ))) ASC  
	)
	--object_id;class_id;class_name;edge_id;geom;confidence;
	SELECT  
		row_number() over() AS qgis_id
		,gid AS object_id 
		,classification_id AS class_id
		, classification AS class_name
		,   edge_id 
		, ST_Force2D(ST_SnapToGrid(ST_SimplifyPreserveTopology(ST_Buffer(ST_Buffer( geom ,2,'quad_segs=4'),-2, 'quad_segs=4'),0.1),0.001)  )  AS geom 
		, ST_Force2D(sline_to_edge) AS sline_to_edge
		,weighted_confidence AS confidence 
	FROM map  ; 

 
---exporting the objects to file
--#object_id;class_id;class_name;edge_id;geom;confidence;
 
COPY ( 
	SELECT regexp_split_to_table(sub.t, E'\\n') 
	FROM (
	 
	( SELECT '##File describing the objects used by snapping' as t
	UNION ALL SELECT '#object_id::int         : unique id of the object'
	UNION ALL SELECT '#class_id::int          : class id of the object'
	UNION ALL SELECT '#class_name::string     : class name of the object'
	UNION ALL SELECT '#edge_id::int           : the objects refers to this edge'
	UNION ALL SELECT '#geom::WKT              : WKT string describing the geom of the object'
	UNION ALL SELECT '#confidence::double     : confidence about this object'
	UNION ALL SELECT '#object_id;class_id;class_name;edge_id;geom;confidence;' 
	)
		UNION ALL 
	( SELECT 	object_id ||';'||
		 class_id||';'||
		 class_name ||';'||
		 edge_id  ||';'||
		  ST_Astext(geom) ||';'||
		 confidence 
	FROM obj_for_output_in_export_area 
	ORDER BY object_id ASC )
	) AS sub
)
TO '/media/sf_E_RemiCura/PROJETS/snapping/data/data_in_reduced_export_area/object_in_export_area.csv';

 