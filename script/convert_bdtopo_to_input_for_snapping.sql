﻿-----------------------------------------------------------
--
--Rémi-C , Thales IGN
--09/2014
--
--This script converts a topological network (edges beign polyline) to a topological network (edge being lines)
--
-----
--details
-----
---NOTE :  
--	doing so we want to conserve full history (what was the number of previous edges/node...)
--	
--warning : uses SFCGAL
-----------------------------------------------------------


CREATE SCHEMA IF NOT EXISTS network_for_snapping; 

SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology, public; 

--creating a table to limit the amount of computing to a given part of the network

/*
DROP TABLE IF EXISTS network_limit ;
CREATE TABLE network_limit(
gid serial PRIMARY KEY,
limit_geom geometry(polygon,932011)) ;
INSERT INTO network_limit (limit_geom) VALUES (ST_GeomFromtext('POLYGON((2065 21785,1458 21241,1838 20694,2403 20590,2597 20765,2594 21393,2406 21757,2065 21785))',932011)); 
*/



--break edges into pair of succesive points
--we give new node_id over 1000000 to nodes resulting from breaking

	DROP TABLE IF EXISTS new_successive_points; 
	CREATE TABLE new_successive_points AS 
	WITH  segmentize AS (
		SELECT  edge_id
			, edge_data.geom --  ST_Segmentize(geom ,5) AS geom
			, start_node, end_node, ST_NumPoints(edge_data.geom)  AS num_points
		 FROM edge_data  , network_limit AS ref_area
		WHERE ST_DWITHIN(geom,limit_geom,40)
	)
	,input_data AS ( --getting all the edges we are going to break into segments .
		--no need to break edge if edge is already a segment and not a polyline!  (hence the tetst on point number)
		SELECT  edge_id, geom , start_node, end_node, ST_NumPoints(geom)  AS num_points
		 FROM segmentize
		 WHERE num_points >2 
		--LIMIT  100 
	)
	 ,points AS ( --we break each edge into successiv points
		SELECT edge_id , dmp.path AS seg_id, dmp.geom ,num_points
		FROM input_data, ST_DumpPoints(geom) AS dmp 
	)
	 ,first_last AS ( --we get the information if the point is the first or the last of the edge. 
				--This is important, because if this is the case, the point will not be a new node but instead the old node that already exists betwen edges in bdtopo_topogical.node
		SELECT * , seg_id[1]= 1 AS is_first, seg_id[1]  = num_points AS is_last 
		FROM points
	) 
	,getting_old_node_id AS ( -- we genereate a new node id : if the point is already a node, keep the correct node_id (start_node or end_node).
			--else, generate a new id starting from 1000000
			-- the variable node_in_intersection keeps tracks of which point was already a node or not
		SELECT nni.edge_id, nni.seg_id[1] As seg_ordinality
			,   COALESCE(ed1.start_node, ed2.end_node, 1000000 + row_number() over()) as new_node_id
			, is_first OR is_last AS node_in_intersection
			,nni.geom
		FROM first_last as nni
			LEFT OUTER JOIN edge_data as ed1 ON (nni.edge_id = ed1.edge_id AND nni.is_first = TRUE)
			LEFT OUTER JOIN edge_data as ed2 ON (nni.edge_id = ed2.edge_id AND nni.is_last = TRUE)
	)
	SELECT *
	FROM getting_old_node_id ;
 

	DROP TABLE IF EXISTS nodes_for_output ; 
	CREATE TABLE nodes_for_output AS  
		--this is the sum of old nodes and newly created nodes. This should be output
			SELECT node_id
				,  ST_X( geom) AS X, ST_Y(geom) AS Y, ST_Z(geom) AS Z
				, 1 AS is_in_intersection
				, geom::geometry(pointZ,932011) AS node_geom
			FROM bdtopo_topological.node , network_limit AS ref_area
			WHERE ST_DWITHIN(geom,limit_geom,40)
		UNION ALL --in theory , we can use UNION ALL. This is a security to enforce no duplicates. We don't care about performance here
			SELECT new_node_id as node_id
				, ST_X( geom) AS X, ST_Y(geom) AS Y, ST_Z(geom) AS Z
				, 0 AS is_in_intersection
				,geom::geometry(pointZ,932011) AS node_geom
			FROM new_successive_points
			WHERE node_in_intersection = FALSE;
			--60 k lines

	CREATE INDEX ON  nodes_for_output USING GIST(ST_SetSRID(ST_MakePoint(X,Y,Z),932011));
	CREATE INDEX ON  nodes_for_output USING GIST(node_geom) ; 
	CREATE INDEX ON nodes_for_output (node_id) ;


	DROP TABLE IF EXISTS edges_for_output ; 
	CREATE TABLE edges_for_output AS  
	WITH successive_points AS (--we get all the successive nodes from getting_old_node_id to create pair of successiv nodes
		SELECT edge_id, seg_ordinality 
			,  new_node_id as start_node
			, lead(new_node_id,1) OVER (PARTITION BY edge_id ORDER BY  seg_ordinality) AS end_node 
		FROM new_successive_points  	
	)
	,new_edges AS (--a pair of nodes is a new edge. We also get the width information.	
			-- we generate one new edge too much per old edge, we have to remove it using the "WHERE sp.end_node IS NOT NULL"
		SELECT 1000000+ row_number() over()  AS edge_id 
			,sp.start_node
			,sp.end_node 
			,r.largeur AS width
			,sp.edge_id AS old_edge_id
		FROM successive_points AS sp
			LEFT OUTER JOIN edge_data AS ed ON (ed.edge_id = sp.edge_id)
			LEFT OUTER JOIN road AS  r ON (ed.ign_id = r.id)
		WHERE sp.end_node IS NOT NULL 
	)
	,edges AS (
		 SELECT ed.edge_id, ed.start_node, ed.end_node, r.largeur AS width, ed.edge_id AS old_edge_id
		FROM edge_data AS ed
		LEFT OUTER JOIN road AS  r ON (ed.ign_id = r.id)
		WHERE ST_NumPoints(ed.geom)<=2
	UNION --should be union all , security against duplicates.
		SELECT  edge_id
				, start_node
				, end_node 
				, width 
				,old_edge_id 
		FROM new_edges
	)
	SELECT edge_id
				, start_node
				, end_node 
				, width 
				,old_edge_id
				, ST_SetSRID(ST_MakeLine(ST_MakePoint(nou1.X,nou1.Y,nou1.Z),ST_MakePoint(nou2.X,nou2.Y,nou2.Z)),932011)::geometry(linestringZ,932011) AS geom
	FROM edges AS eou
		LEFT OUTER JOIN nodes_for_output AS nou1 ON (eou.start_node = nou1.node_id)
		LEFT OUTER JOIN nodes_for_output AS nou2 ON (eou.end_node = nou2.node_id) ;
		
	CREATE INDEX ON edges_for_output USING GIST(geom) ;
	 CREATE INDEX ON edges_for_output (edge_id) ;
	 CREATE INDEX ON edges_for_output (start_node) ;
	  CREATE INDEX ON edges_for_output (end_node) ;
 

--now we restrain ourselve to an area for test :
	--this is delimitated by the table def_zone_export

/*
DROP TABLE IF EXISTS def_zone_export; 
	CREATE TABLE def_zone_export 
	( 	gid SERIAL PRIMARY KEY
		,id bigint
		,geom geometry(POLYGON,931008)
		, center_in_RGF93 geometry(POINT,932011) );   
	INSERT INTO def_zone_export (id, geom) 
		VALUES( 1 , ST_GeomFromText('POLYGON((650907.6  6860870.6 ,650956.9  6860895.8 ,651036.0  6860757.1 ,650983.7  6860736.6 ,650907.6  6860870.6 ))',931008) );

 	INSERT INTO def_zone_export (id, geom)  --all_of_acquisition
 		VALUES( 1 , ST_GeomFromText('POLYGON((650637.4 6861097.2,650909.8 6861534.1,651068.6 6861638.3,651365.1 6861697.3,651410.3 6861270.3,651544.3 6861218.2,651498.6 6861178.5,651288.1 6861198.7,650863.3 6861051.6,651004.5 6860830.8,651063.1 6860739.1,651421.5 6860706.1,651347.3 6860611.3,651030.6 6860666.2,650657 6861002.8,650637.4 6861097.2))',931008) );
	UPDATE def_zone_export SET center_in_RGF93 = ST_Centroid(ST_Transform(geom,932011)) ; 
*/

		--importing the sidewalk observation (weighted points ) here
/*
	DROP MATERIALIZED VIEW IF EXISTS weighted_sidewalk_observation_point ;
	CREATE MATERIALIZED VIEW weighted_sidewalk_observation_point AS 
	SELECT *
	FROM public.dblink( 
		 'hostaddr=127.0.0.1 port=5433 dbname=TerraMobilita user=postgres password=postgres'::text
		, 'SELECT * FROM cmm.lines_to_weighted_points '::text
		,TRUE)  AS   f(
			 qgis_id bigint,
			  area_id bigint,
			  seg_id integer,
			  geom geometry,
			  weight numeric
			 ) ; 
	CREATE INDEX ON weighted_sidewalk_observation_point (qgis_id) ; 
	CREATE INDEX ON weighted_sidewalk_observation_point (area_id) ; 
	CREATE INDEX ON weighted_sidewalk_observation_point (seg_id) ; 
	CREATE INDEX ON weighted_sidewalk_observation_point (weight) ; 
	CREATE INDEX ON weighted_sidewalk_observation_point USING GIST  (ST_Transform(geom,932011) ); 



	DROP TABLE IF EXISTS weighted_sidewalk_observation_line ;
	CREATE TABLE weighted_sidewalk_observation_line AS 
	SELECT *
	FROM public.dblink( 
		 'hostaddr=127.0.0.1 port=5433 dbname=TerraMobilita user=postgres password=postgres'::text
		, '	SELECT row_number() over() as gid, area_id,  dmp.geom AS geom   
			FROM cmm.filtered_lines_final, ST_Dump(ST_Transform(area_skeleton,932011) ) AS dmp
			'::text
		,TRUE)  AS   f( 
			observation_id bigint, 
			  area_id bigint,
			  geom geometry(linestring, 932011)  
			 ) ; 
	ALTER TABLE weighted_sidewalk_observation_line ADD PRIMARY KEY (observation_id ); 
	CREATE INDEX ON weighted_sidewalk_observation_line USING GIST  (geom ); 
  
*/
	-------------------------
	DROP TABLE IF EXISTS trottoir_cut_into_pieces ;  
	CREATE TABLE trottoir_cut_into_pieces AS 
	WITH trottoir AS ( --filtering to keep only ground cornerstones
		SELECT t.gid, t.geom 
		FROM odparis_corrected.trottoir as t  , def_zone_export as dfz 
		WHERE niveau = 'Surface'
			AND libelle = 'Bordure'  
			AND ST_Area( Box2D(t.geom) ) >20
			AND ST_Length( t.geom) >14
			AND ST_INtersects( ST_Transform(t.geom,932011), ST_Transform(dfz.geom,932011))=TRUE 
	)
	SELECT row_number() over() as qgis_id, 1.0 as weight, 1.0 as confidence, gid, dmp.geom::geometry(point,932011)  as geom
	FROM trottoir,ST_DumpPoints(ST_Transform(ST_Segmentize(ST_SimplifyPreserveTopology(ST_Segmentize( geom, 4.0),0.5), 4.0),932011)) AS dmp ; 

	CREATE INDEX ON trottoir_cut_into_pieces USING GIST (geom) ;
	-- CREATE INDEX ON  trottoir_cut_into_pieces USING GIST (ST_Transform( geom,932011)) ; 
	CREATE INDEX ON trottoir_cut_into_pieces (gid); 
	--qgis_id,geom, weight ;
	 




	DROP TABLE IF EXISTS edges_for_output_in_export_area CASCADE; 
	CREATE TABLE edges_for_output_in_export_area AS 
		--we keep only the edges whose both ends are exported nodes.
	SELECT DISTINCT edge_id
				, start_node
				, end_node 
				, width 
				, old_edge_id 
				, eou.geom::geometry(linestringZ,932011)
	FROM edges_for_output as eou ,  def_zone_export as dfz
	WHERE ST_Intersects(eou.geom, ST_Transform((dfz.geom),932011)  ) 
		AND edge_id is not null;
	ALTER TABLE edges_for_output_in_export_area ADD PRIMARY KEY (edge_id) ;  
	CREATE INDEX ON edges_for_output_in_export_area (start_node);
	CREATE INDEX ON edges_for_output_in_export_area (end_node);
	CREATE INDEX ON edges_for_output_in_export_area USING GIST(geom);
  
	
	DROP TABLE IF EXISTS nodes_for_output_in_export_area CASCADE; 
	CREATE TABLE nodes_for_output_in_export_area AS 
		--we keep only the nodes inside the export area
	WITH node_id AS (SELECT start_node  AS node_id
	FROM edges_for_output_in_export_area
	UNION
	SELECT end_node  AS node_id
	FROM edges_for_output_in_export_area
	)
	SELECT nfo.node_id
		, nfo.X AS X
		,nfo.Y AS Y 
		,COALESCE(nfo.Z,0) AS Z
		,nfo.is_in_intersection 
		, ST_SetSRID(ST_MakePoint(X,Y,Z),932011)::geometry(POINTZ,932011) as node_geom
	FROM node_id
		NATURAL JOIN nodes_for_output AS nfo;

/*
	-- small trigger on node for output so XYZ are always in sync
	DROP FUNCTION IF EXISTS rc_update_XYZ_when_geom_change(  ) CASCADE; 
	CREATE OR REPLACE FUNCTION rc_update_XYZ_when_geom_change(  ) RETURNS  trigger  AS $BODY$  
	DECLARE  
		BEGIN 
			IF geometry_eq(NEW.node_geom, OLD.node_geom) THEN NEW.node_geom :=ST_SetSRID(ST_MakePoint(NEW.X,NEW.Y,NEW.Z),932011) ; 
			ELSE
			NEW.X := ST_X(node_geom) ;  NEW.Y := ST_Y(node_geom) ;  NEW.Z := ST_Z(node_geom) ; END IF;
			RETURN NEW;
		END ;
	$BODY$ LANGUAGE plpgsql VOLATILE;
	DROP VIEW IF EXISTS nodes_for_output_in_export_area_qgis ; 
	CREATE VIEW   nodes_for_output_in_export_area_qgis AS
	SELECT node_id, X,Y,Z,is_in_intersection, ST_Force2D(node_geom)::geometry(point,932011) as node_geom FROM nodes_for_output_in_export_area; 

	DROP TRIGGER IF EXISTS  rc_update_XYZ_when_geom_change ON nodes_for_output_in_export_area; 
	CREATE  TRIGGER rc_update_XYZ_when_geom_change   BEFORE  UPDATE   ON nodes_for_output_in_export_area FOR EACH ROW  
	EXECUTE PROCEDURE rc_update_XYZ_when_geom_change(); 
*/

	-- creating a tbale for user override 
	-- DROP TABLE IF EXISTS sidewalk_override ; 
	CREATE TABLE IF NOT EXISTS sidewalk_override  (
	gid serial primary key,
	user_point geometry(POINT,932011) ,
	weight float DEFAULT 10
	);
	--TRUNCATE sidewalk_override ; 
	-- INSERT INTO sidewalk_override (user_point) VALUES (ST_GeomFromText('POINT(1950 2140)',932011)) ; 
 
	-- CREATE INDEX ON sidewalk_override USING GIST(user_point) ; 
	--UPDATE sidewalk_override SET weight = 50; 

	--getting the observations inside the area 
  
	DROP TABLE IF EXISTS obs_for_output_in_export_area; 
	CREATE TABLE obs_for_output_in_export_area AS  
	WITH map AS (--for each edge, we get observation closer than 5 meters 
		(SELECT DISTINCT ON (oia.qgis_id) qgis_id, area_id, seg_id, oia.geom, weight, eg.edge_id
		FROM def_zone_export as dfz , weighted_sidewalk_observation_point  as oia , edges_for_output_in_export_area AS eg
		WHERE ST_DWithin( ST_Transform(oia.geom,932011), eg.geom,4+width)=TRUE
			AND ST_DWithin( ST_Transform(oia.geom,932011), eg.geom,20)=TRUE
			AND ST_WITHIN( ST_Transform(oia.geom,932011),  ST_Transform(dfz.geom,932011) ) = TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface,ST_Transform(oia.geom,932011),1)=TRUE
			)
		ORDER BY oia.qgis_id ASC, ST_Distance(ST_Transform(oia.geom,932011),eg.geom ) ASC )

		UNION ALL 

		(SELECT DISTINCT ON (oia.gid) oia.gid + (SELECT max(qgis_id) FROM weighted_sidewalk_observation_point),NULL, NULL, user_point AS geom , weight, eg.edge_id
		FROM def_zone_export as dfz , sidewalk_override  as oia , edges_for_output_in_export_area AS eg
		WHERE ST_DWithin(  oia.user_point , eg.geom,4+width)=TRUE
			AND ST_DWithin( oia.user_point , eg.geom,20)=TRUE
			AND ST_WITHIN(  oia.user_point ,  ST_Transform(dfz.geom,932011) ) = TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface, oia.user_point ,1)=TRUE
			)
		ORDER BY oia.gid ASC, ST_Distance( oia.user_point ,eg.geom ) ASC) 

		UNION ALL
		(SELECT DISTINCT ON (oia.qgis_id) oia.qgis_id + (SELECT max(qgis_id) FROM weighted_sidewalk_observation_point) +  (SELECT count(*) FROM sidewalk_override) 
			,NULL, NULL,  oia.geom , weight, eg.edge_id
		FROM def_zone_export as dfz , trottoir_cut_into_pieces  as oia , edges_for_output_in_export_area AS eg
		WHERE ST_DWithin(  oia.geom , eg.geom,4+width)=TRUE
			AND ST_DWithin( oia.geom , eg.geom,20)=TRUE
			AND ST_WITHIN(  oia.geom ,  ST_Transform(dfz.geom,932011) ) = TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface, oia.geom ,3)=TRUE
			)
		ORDER BY oia.qgis_id ASC, ST_Distance( oia.geom ,eg.geom ) ASC) 
	)
	SELECT row_number() over() AS obs_id 
		, edge_id
		, X,Y,Z
		, 1::float AS confidence, COALESCE(weight,0) AS weight
		, ST_SetSRID(ST_MakePoint(X,Y), 932011)::geometry(point,932011) AS geom
	FROM map, ST_Transform( map.geom,932011) As tgeom
		 , ST_X(tgeom) AS X
			, ST_Y(tgeom) AS Y
			, COALESCE(ST_Z(tgeom),0) AS Z ;
 
	ALTER TABLE obs_for_output_in_export_area ADD  PRIMARY KEY (obs_id) ; 

/*
 	COPY  nodes_for_output_in_export_area TO '/media/sf_USB_storage/PROJETS/snapping/data/data_in_reduced_export_area/nodes_for_output_in_export_area.csv'	 (FORMAT CSV, DELIMITER   ';') ;
 	 COPY( SELECT edge_id, start_node, end_node , width,old_edge_id FROM  edges_for_output_in_export_area )TO '/media/sf_USB_storage/PROJETS/snapping/data/data_in_reduced_export_area/edges_for_output_in_export_area.csv'	 (FORMAT CSV, DELIMITER   ';') ;
	COPY ( SELECT  obs_id , edge_id , X , Y , Z  , confidence, weightobs_for_output_in_export_area) 
	TO '/media/sf_USB_storage/PROJETS/snapping/data/data_in_reduced_export_area/obs_for_output_in_export_area.csv' (FORMAT CSV, DELIMITER   ';') ;
*/

	
--visualization
	--we want to visualize edges, nodes and observations alongs with pairings

	--visaulize node
	DROP TABLE IF EXISTS nodes_for_visu_in_export_area; 
	CREATE TABLE nodes_for_visu_in_export_area AS 
	SELECT node_id, ST_SetSRID(ST_MakePoint(X::float,Y::float,Z::float),932011) AS geom, is_in_intersection
	FROM nodes_for_output_in_export_area;  

	CREATE INDEX ON nodes_for_visu_in_export_area USING GIST(geom);
	CREATE INDEX ON nodes_for_visu_in_export_area (node_id);

	--visualize edges :
	DROP TABLE IF EXISTS edges_for_visu_in_export_area; 
	CREATE TABLE edges_for_visu_in_export_area AS 
	SELECT edge_id, ST_MakeLine(n1.geom,n2.geom) AS geom, width
	FROM edges_for_output_in_export_area AS ed
		LEFT OUTER JOIN nodes_for_visu_in_export_area AS n1 ON (n1.node_id = ed.start_node)
		LEFT OUTER JOIN nodes_for_visu_in_export_area AS n2 ON (n2.node_id = ed.end_node);

	--visualize observations :
	DROP TABLE IF EXISTS obs_for_visu_in_export_area; 
	CREATE TABLE obs_for_visu_in_export_area AS 
	SELECT obs_id, obf.edge_id, sgeom AS geom
		, confidence, weight
		, ST_SHortestLine(sgeom, ed.geom) AS sline_to_edge
	FROM obs_for_output_in_export_area AS obf
		LEFT JOIN edges_for_visu_in_export_area AS ed USING(edge_id)
		,ST_SetSRID(ST_MakePoint(X::float,Y::float,Z::float),932011) AS sgeom;  
 

	--visualize observation relating to street gen
	DROP TABLE IF EXISTS obs_r_street_gen_in_export_area; 
	CREATE TABLE obs_r_street_gen_in_export_area AS 
	SELECT DISTINCT ON (obs_id) obs_id,  sgeom AS geom
		, confidence, weight
		, ST_SHortestLine(sgeom, ST_ExteriorRing(ra.section2_surface)) AS sline_to_edge
		,ST_Intersects(ra.section2_surface, sgeom) as is_inside
	FROM obs_for_output_in_export_area AS obf
		,ST_SetSRID(ST_MakePoint(X::float,Y::float,Z::float),932011) AS sgeom 
		INNER JOIN  street_amp.result_axis as ra ON ST_DWithin(ra.section2_surface, sgeom, 5)   
	ORDER BY obs_id , ST_Distance(ra.section2_surface, sgeom) ASC ;



-- computing target slope information
 
	
	DROP TABLE IF EXISTS slope_for_output_in_export_area CASCADE; 
	CREATE TABLE slope_for_output_in_export_area AS  
	WITH map AS (  
		SELECT DISTINCT ON (oia.observation_id, dmp.path) 
			 observation_id, dmp.geom AS seg_geom
			, dmp.path as seg_order
			, eg.edge_id
		
		FROM def_zone_export as dfz , weighted_sidewalk_observation_line   AS oia  , 
			edges_for_output_in_export_area AS eg
			--,edge_geom AS eg
			,rc_lib.rc_DumpSegments(ST_Segmentize(oia.geom, 2.0)) as dmp  
		WHERE ST_DWithin(oia.geom, ST_transform(dfz.geom ,932011),30 ) = TRUE
			AND ST_DWithin(oia.geom,eg.geom,30 ) = TRUE
			AND ST_DWithin( dmp.geom, eg.geom,30)=TRUE
			AND ST_DWithin( dmp.geom, eg.geom,4+width)=TRUE 
			AND ST_WITHIN(dmp.geom,  ST_transform(dfz.geom ,932011)) = TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface,dmp.geom,3)=TRUE
			)
		ORDER BY oia.observation_id, dmp.path, ST_Distance(dmp.geom,eg.geom ) ASC 
	)
	, angles AS (
		SELECT map.*
			, CASE WHEN degrees >180-20 THEN degrees  - 180 ELSE degrees END AS angle
			, degrees
			, ST_Length(seg_geom) as weight
		FROM map
			, mod(CAST(degrees(ST_Azimuth(ST_StartPoint(seg_geom), ST_EndPoint(seg_geom)))AS INT),180) as degrees
	)
	, weighted_median AS (
		SELECT edge_id
			, rc_py_weighted_median(
				array_agg(angle order by observation_id, seg_order )
				,array_agg(weight order by observation_id, seg_order ) 
			)::int as weighted_median
		FROM angles
		GROUP BY edge_id
	) 
	, angle_and_median AS (
		SELECT edge_id, angle, weight, abs(angle - weighted_median) < 20 AS to_be_used
			-- ,  sum(angle* weight) / (SELECT sum(weight) FROM angles) AS target_slope 
			-- , diff, mod
		FROM angles
			LEFT OUTER JOIN weighted_median USING (edge_id)  
		-- WHERE --  abs(angle - weighted_median) < 20 = TRUE AND
		-- edge_id = 1000080
	)
	, filtered AS (
		SELECT edge_id, sum(angle*weight) / sum(weight) AS target_slope
			, variance(angle) AS variance
			, sum(weight) AS weights
		FROM angle_and_median
		WHERE to_be_used = TRUE
		GROUP BY edge_id 
	)	 
	SELECT edge_id,mod(target_slope::int+180,180) AS target_slope, 1 - sqrt(variance) / (180 ) as confidence, weights AS weight
	FROM filtered 
	WHERE variance IS NOT NULL  ;-- remove case with only one observation;  
	
	ALTER TABLE slope_for_output_in_export_area ADD CONSTRAINT distfk FOREIGN KEY (edge_id) REFERENCES edges_for_output_in_export_area (edge_id) MATCH FULL;

	DROP VIEW IF EXISTS slope_for_output_in_export_area_v ; 
	CREATE VIEW slope_for_output_in_export_area_v AS
	SELECT edge_id, ST_Centroid(eg.geom)::geometry(point,932011) AS geom, target_slope, confidence, weight
	FROM slope_for_output_in_export_area AS fa
		LEFT OUTER JOIN edges_for_output_in_export_area AS eg USING(edge_id)  ; 


		 

---test to export the whole file complete, no 3 separated files :

 
 
COPY (
	 SELECT '#num nodes;num edge;num observation'
		UNION ALL 
	 SELECT (SELECT count(*) FROM nodes_for_output_in_export_area) || ';' 
		||  (SELECT count(*) FROM edges_for_output_in_export_area) || ';' 
		|| (SELECT count(*) FROM obs_for_output_in_export_area) || ';'
		|| (SELECT count(*) FROM slope_for_output_in_export_area)
		UNION ALL 
	 SELECT '#node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int'
		UNION ALL 
	SELECT node_id ||';'||x||';'||y||';'||z||';'||z||';'||is_in_intersection||';' 
	FROM nodes_for_output_in_export_area
		UNION ALL 
	 SELECT '#edge_id::int;start_node::int;end_node::int;width::double'
		UNION ALL 
	SELECT edge_id ||';'||start_node||';'||end_node||';'||width
	FROM edges_for_output_in_export_area
		UNION ALL 
	 SELECT '#obs_id::int;edge_id::int;X::double;Y::double;Z::double;confidence::double;weight::double'
		UNION ALL 
	SELECT obs_id ||';'||edge_id||';'||x||';'||y||';'||z||';'||confidence||';'||weight
	FROM obs_for_output_in_export_area
		UNION ALL 
	SELECT '#edge_id::int;slope::double;confidence::double;weight::double'
		UNION ALL 
	SELECT edge_id||';'|| target_slope||';'|| confidence||';'||weight
	FROM slope_for_output_in_export_area
)
TO '/media/sf_USB_storage/PROJETS/snapping/data/data_in_reduced_export_area/full_area.csv' ;


SELECT degrees(ST_Azimuth(st_startpoint(geom), st_endpoint(geom)))
FROM edges_for_output_in_export_area
WHERE edge_id = 1000041 ; 

SELECT degrees(-0.15708)

DROP TABLE IF EXISTS tem ;
CREATE TABLE tem (
gid serial PRIMARY KEY,
geom geometry(point,932011)) ;  

TRUNCATE tem ;
INSERT INTO tem (geom) VALUES  
	(ST_GeomFromtext('POINT(1957.73  21288)',932011))
	, (ST_GeomFromtext('POINT(1961.52 21342.1)',932011)) ; 
 


 