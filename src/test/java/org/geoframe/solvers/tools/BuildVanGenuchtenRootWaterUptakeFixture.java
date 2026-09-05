package org.geoframe.solvers.tools;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler;
import org.hortonmachine.dbs.compat.ADb;
import org.hortonmachine.dbs.compat.EDb;
import org.hortonmachine.dbs.utils.SqlName;

/**
 * One-off generator for {@code RichardsCoupled_VG_RootUptake.gpkg}, used by
 * {@code TestVanGenuchtenRootWaterUptakeGpkg}. Copies the real, already-committed
 * {@code RichardsCoupled_VG.gpkg} fixture (same 3-layer soil, same real rainfall
 * BC {@code TestVanGenuchtenGpkg} already uses - no CSV source exists for it,
 * only the compiled binary) and adds {@link Whetgeo1DInputsHandler#COL_SWRC_THETA_WP}/
 * {@link Whetgeo1DInputsHandler#COL_SWRC_THETA_FC} columns to its
 * {@code swrc_parameters} table.
 *
 * <p>
 * Unlike the rest of this fixture (real soil/rainfall data), the thetaWP/thetaFC
 * values here are synthetic placeholders, not real site data - there's no
 * historical scenario to match, this is a new solver-mechanics test, not a port.
 * Placed at a generic, physically sane fraction of each layer's own
 * {@code thetaR..thetaS} range (thetaWP at 30%, thetaFC at 70%) so the clamp has
 * real room to matter without being trivially always-on or never-on.
 */
public class BuildVanGenuchtenRootWaterUptakeFixture {

	public static void main(String[] args) throws Exception {
		String srcPath = "src/test/resources/input/gpkg/RichardsCoupled_VG.gpkg";
		String outPath = "src/test/resources/input/gpkg/RichardsCoupled_VG_RootUptake.gpkg";

		Files.copy(new File(srcPath).toPath(), new File(outPath).toPath(), StandardCopyOption.REPLACE_EXISTING);

		try (ADb db = EDb.GEOPACKAGE.getDb()) {
			db.open(outPath);

			SqlName table = SqlName.m(Whetgeo1DInputsHandler.TABLE_SWRC_PARAMETERS);
			db.executeInsertUpdateDeleteSql(
					"ALTER TABLE " + table.fixedDoubleName + " ADD COLUMN " + Whetgeo1DInputsHandler.COL_SWRC_THETA_WP
							+ " REAL");
			db.executeInsertUpdateDeleteSql(
					"ALTER TABLE " + table.fixedDoubleName + " ADD COLUMN " + Whetgeo1DInputsHandler.COL_SWRC_THETA_FC
							+ " REAL");

			String updateSql = "UPDATE " + table.fixedDoubleName + " SET "
					+ Whetgeo1DInputsHandler.COL_SWRC_THETA_WP + " = "
					+ Whetgeo1DInputsHandler.COL_SWRC_THETAR + " + 0.3 * ("
					+ Whetgeo1DInputsHandler.COL_SWRC_THETAS + " - " + Whetgeo1DInputsHandler.COL_SWRC_THETAR + "), "
					+ Whetgeo1DInputsHandler.COL_SWRC_THETA_FC + " = "
					+ Whetgeo1DInputsHandler.COL_SWRC_THETAR + " + 0.7 * ("
					+ Whetgeo1DInputsHandler.COL_SWRC_THETAS + " - " + Whetgeo1DInputsHandler.COL_SWRC_THETAR + ")";
			db.executeInsertUpdateDeleteSql(updateSql);
		}
		System.out.println("Wrote " + outPath);
	}
}
