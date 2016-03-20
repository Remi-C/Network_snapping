
/*************************************
**Remi C, 2016, Thales IGN
**
*****************************************/

/*slight addon: how to find target slope for an edge*/



SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology, public; 


DROP FUNCTION IF EXISTS rc_py_weighted_median( val FLOAT[] , weight FLOAT[] );
CREATE FUNCTION rc_py_cluster_points_dbscan( val FLOAT[] , weight FLOAT[], OUT weighted_median FLOAT)
AS $$"""
This function take a list of values and wieght, are return the weighted median
require weightedstats
"""

import numpy as np
import weightedstats as ws
val_ = np.array( val, dtype=np.float)
weight_ = np.array( weight, dtype=np.float) 

weighted_median = ws.numpy_weighted_median(val_, weights=weight_)
 
return weighted_median
$$ LANGUAGE plpythonu IMMUTABLE STRICT; 

WITH all_observations AS (

)

WITH observations AS(
	SELECT
		edge_id, 
		,observation_id as polyline_id, dmp.path
		, dmp.geom as line_segment
		, confidence
		, st_length(geom) as weight
		, degrees(ST_Azimuth(ST_StartPoint(geom), ST_EndPoint(geom)))::int % 180 as angle
	FROM observations, rc_lib.rc_DumpSegments(lines) as dmp
)
, weighted_median AS (
	SELECT rc_py_weighted_median(array_agg(angle order by ), )
	FROM observations
)


------------------------------
DROP TABLE IF EXISTS slope_for_output_in_export_area; 
	CREATE TABLE slope_for_output_in_export_area AS  
	WITH mapping AS (  
		SELECT DISTINCT ON (oia.observation_id, dmp.path) 
			oia.observation_id, seg_geom
			, dmp.path as seg_order
			, eg.edge_id
		
		FROM observation   AS oia  , def_zone_export as dfz ,edge_geom AS eg
			,rc_lib.rc_DumpSegments(lines) as dmp, ST_Transform(dmp.geom,932011) AS seg_geom
		WHERE ST_DWithin( seg_geom, eg.edge_geom,30)=TRUE
			AND ST_DWithin( seg_geom, eg.edge_geom,4+width)=TRUE 
			AND ST_WITHIN(oia.geom,  dfz.geom ) = TRUE 
			AND NOT EXISTS (
				SELECT 1
				FROM street_amp.result_intersection as ri
				WHERE ST_DWithin(ri.intersection_surface,seg_geom,3)=TRUE
			)
		ORDER BY oia.observation_id, dmp.path, ST_Distance(seg_geom,eg.edge_geom ) ASC
	)
	, angles AS (
		SELECT *
			, degrees(ST_Azimuth(ST_StartPoint(seg_geom), ST_EndPoint(seg_geom)))::int % 180 as angle
			, ST_Length(seg_geom) as weight
		FROM mapping
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
	, filtered_angle AS(
		SELECT edge_id, sum(angle* weight) / (SELECT sum(weight) FROM angles) AS target_slope
			, variance(mod_angled) as variance
			, sum(weight) AS weights
		FROM weighted_median, mod(weighted_median + 90 , 180) AS mod_median
			,angles, mod(angle + weighted_median + 90 , 180) AS mod_angled
		WHERE mod_angled BETWEEN mod_median - 20, mod_median + 20 
	)
	SELECT edge_id, target_slope, 1 - variance / (20 + 20 ) as confidence, weights AS weight
	FROM filtered_angle ; 
 



